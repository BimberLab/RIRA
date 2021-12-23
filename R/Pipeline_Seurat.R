library(mlr3verse)



#' @title Creates a binary classifier to classify cells within a Seurat object
#'
#' @description Creates a binary classifier to classify cells via TrainAllModels to iterate over all cell types.
#' @param training_matrix A counts or data slot provided by TrainAllModels
#' @param celltype The celltype (provided by TrainAllModels) used as classifier's positive prediction
#' @param hyperparameter_tuning logical that determines whether or not hyperparameter tuning should be performed.
#' @param learner The mlr3 learner that should be used. Currently fixed to "classif.ranger" if hyperparameter tuning is FALSE. Otherwise, "classif.xgboost" and "classif.ranger" are supported.
#' @param inner_resampling The resampling strategy that is used for hyperparameter optimization. Holdout ("hout" or "holdout") and cross validation ("cv" or "cross-validation") are supported.
#' @param outer_resampling The resampling strategy that is used to determine overfitting. Holdout ("hout" or "holdout") and cross validation ("cv" or "cross-validation") are supported.
#' @param inner_folds The number of folds to be used for inner_resampling if cross-valdiation is performed.
#' @param outer_folds The number of folds to be used for outer_resampling if cross-valdiation is performed.
#' @param inner_ratio The ratio of training to testing data to be used for inner_resampling if holdout resampling is performed.
#' @param outer_ratio The ratio of training to testing data to be used for inner_resampling if holdout resampling is performed.
#' @param n_models The number of models to be trained during hyperparameter tuning. The model with the highest accuracy will be selected and returned.
#' @param parallel Logical determining if futures should be used or not for parallelization
#' @export

TrainModel <- function(training_matrix, celltype, hyperparameter_tuning = F, learner = "classif.ranger", inner_resampling = "cv", outer_resampling = "cv", inner_folds = 4, inner_ratio = 0.8,  outer_folds = 3, outer_ratio = 0.8, n_models = 20, parallel = F){
  GetSeed()
  if (parallel){
    future::plan("multiprocess")
  }
  #fix gene names to conform with Seurat Object processing
  colnames(training_matrix) <- gsub("-", ".", names(training_matrix))
  #create classification matrix 
  classification.data <- training_matrix
  #Trim the celltype_column containing all of the ground truth celltypes 
  #Note, this column is named "celltype" and is not the same as the column named by the passed celltype argument, which changes according to the celltype being predicted
  classification.data <- subset(classification.data, select = -celltype)
  #trim the binarized label/celltype/truth column (should be the last column)
  classification.data <- classification.data[,1:(ncol(classification.data)-1)]
  
  celltype_binary <- training_matrix[,ncol(training_matrix)]
  classification.data[,"celltype_binary"] <- as.factor(celltype_binary)
  colnames(classification.data) <- make.names(colnames(classification.data),unique = T)
  
  task <- TaskClassif$new(classification.data, id = "CellTypeBinaryClassifier", target = "celltype_binary")
  
  if (!hyperparameter_tuning){
    #Use a ranger random forest tree with default parameters and holdout (80% train, 20% test) resampling
    #classification.data has column number = number_of_genes + 1, so the mtry argument uses the full gene matrix
    learner <- lrn("classif.ranger", importance = "permutation", num.trees=500, mtry=(ncol(classification.data)-1), predict_type = "prob")
    set_threads(learner)
    train_set = sample(task$nrow, 0.8 * task$nrow)
    test_set = setdiff(seq_len(task$nrow), train_set)
    model <- learner$train(task, row_ids = train_set)$predict(task, row_ids=test_set)
    confusion <- caret::confusionMatrix(factor(model$response), factor(model$truth))
    return(list(model= learner, metrics=confusion))
    
  } else if (hyperparameter_tuning){
    #Set model-independent values for the autotuner
    measure <- msr("classif.ce")
    terminator <- trm("evals", n_evals = n_models)
    #Define a tuning space 25% as large as the number of models 
    #In the case of sensitive hyperparameters, resolution = 5 allows for a low/medium-low/medium/medium-high/high type parameter space
    tuner <- tnr("grid_search", resolution = 5)
    
    #Define resampling method used for hyperparameter tuning
    if (inner_resampling == "cv" | inner_resampling == "cross-validation"){
      inner_resample <- rsmp("cv", folds = inner_folds)
    } else if (inner_resampling == "hout" | inner_resampling == "holdout"){
      inner_resample <- rsmp("holdout", ratio = inner_ratio)
    } else {
      print("Unknown inner_resampling method provided. Please select one of cross-validation, cv, hout, or holdout")
    }
    #Define resampling method used to determine overfitting
    if (outer_resampling == "cv" | outer_resampling == "cross-validation"){
      outer_resample <- rsmp("cv", folds = outer_folds)
    } else if (outer_resampling == "hout" | outer_resampling == "holdout"){
      outer_resample <- rsmp("holdout", ratio = outer_ratio)
    } else {
      print("Unknown outer_resampling method provided. Please select one of cross-validation, cv, hout, or holdout")
    }
    
    #Set learner and define parameter space
    if(learner == "classif.ranger"){
      #Define learner
      learner <- lrn("classif.ranger", importance = "permutation", predict_type = "prob")
      #Define Ranger Hyperparameter Space (RandomBotv2)
      tune_ps = ps(
        num.trees = p_int(lower = 10, upper = 2000),
        sample.fraction = p_dbl(lower = 0.1, upper = 1),
        respect.unordered.factors = p_fct(levels = c("ignore", "order", "partition")), 
        min.node.size = p_int(lower = 1, upper = 100), 
        splitrule = p_fct(levels = c("gini", "extratrees")),
        num.random.splits = p_int(lower = 1, upper = 100, depends = splitrule == "extratrees")
        )
    } else if (learner == "classif.xgboost"){
      #Update task
      task <- TaskClassif$new(classification.data, id = "CellTypeBinaryClassifier", target = "celltype_binary")
      #Define learner
      learner <- lrn("classif.xgboost", predict_type = "prob")
      #Define XGBoost model's Hyperparameter Space (RandomBotv2)
      tune_ps = ps(
        booster = p_fct(levels = c("gblinear", "gbtree", "dart")),
        nrounds = p_int(lower = 2, upper = 8, trafo = function(x) as.integer(round(exp(x)))),
        eta = p_dbl(lower = -4, upper = 0, trafo = function(x) 10^x), 
        gamma = p_dbl(lower = -5, upper = 1, trafo = function(x) 10^x), 
        lambda = p_dbl(lower = -4, upper = 3, trafo = function(x) 10^x),
        alpha = p_dbl(lower = -4, upper = 3, trafo = function(x) 10^x), 
        subsample = p_dbl(lower = 0.1, upper = 1),
        max_depth = p_int(lower = 1, upper = 15),
        min_child_weight = p_dbl(lower = -1, upper = 0, trafo = function(x) 10^x),
        colsample_bytree = p_dbl(lower = 0.1, upper = 1),
        colsample_bylevel = p_dbl(lower = 0.1, upper = 1),
        rate_drop = p_int(lower = 0, upper = 1, depends = booster == 'dart'),
        skip_drop = p_int(lower = 0, upper = 1, depends = booster == 'dart')
      )
      
      }
    }
    
    #Define the autotuner using the parameter spaces and conditions defined above
    at <- AutoTuner$new(
      learner = learner,
      resampling = inner_resample,
      measure = measure,
      search_space = tune_ps,
      terminator = terminator,
      tuner = tuner
    )
    
    #Train the initial model to optimize hyperparameters
    at$train(task)
    
    #instantiate outer resampling
    outer_resample$instantiate(task)
    #score model using outer resampling strategy & output full resampling to be mined for metrics
    full_resampling <- resample(task, at, outer_resample, store_models = TRUE)
    
    return(list(model = at, metrics = full_resampling))
  }


#' @title Wrapper function for TrainModel to train a suite of binary classifiers for each cell type present in the data
#'
#' @description Wrapper function for TrainModel to train a suite of binary classifiers for each cell type present in the data
#' @param seuratObj The Seurat Object to be updated
#' @param celltype_column The metadata column containing the celltypes. One classifier will be created for each celltype present in this column.
#' @param assay SeuratObj assay containing the desired count matrix/metadata
#' @param slot Slot containing the count data. Should be restricted to counts, data, or scale.data. 
#' @param output_dir The directory in which models, metrics, and training data will be saved. 
#' @param hyperparameter_tuning Logical that determines whether or not hyperparameter tuning should be performed.
#' @param learner The mlr3 learner that should be used. Currently fixed to "classif.ranger" if hyperparameter tuning is FALSE. Otherwise, "classif.xgboost" and "classif.ranger" are supported.
#' @param inner_resampling The resampling strategy that is used for hyperparameter optimization. Holdout ("hout" or "holdout") and cross validation ("cv" or "cross-validation") are supported.
#' @param outer_resampling The resampling strategy that is used to determine overfitting. Holdout ("hout" or "holdout") and cross validation ("cv" or "cross-validation") are supported.
#' @param inner_folds The number of folds to be used for inner_resampling if cross-valdiation is performed.
#' @param outer_folds The number of folds to be used for outer_resampling if cross-valdiation is performed.
#' @param inner_ratio The ratio of training to testing data to be used for inner_resampling if holdout resampling is performed.
#' @param outer_ratio The ratio of training to testing data to be used for inner_resampling if holdout resampling is performed.
#' @param n_models The number of models to be trained during hyperparameter tuning. The model with the highest accuracy will be selected and returned.
#' @param parallel Logical determining if futures should be used or not for parallelization
#' @param verbose Whether or not to print the metrics data for each model after training. 
#' @export


TrainAllModels <- function(seuratObj, celltype_column, assay = "RNA", slot = 'data', output_dir = "./classifiers", hyperparameter_tuning = F, learner = "classif.ranger", inner_resampling = "cv", outer_resampling = "cv", inner_folds = 4, inner_ratio = 0.8,  outer_folds = 3, outer_ratio = 0.8, n_models = 20, parallel = F, verbose = T){
  GetSeed()
  if (endsWith(output_dir, "/")){
    print("Please remove the trailing slash in output_dir")
    stop()
  }
  
  #Read the raw data from a seurat object and parse into an mlr3-compatible labeled matrix
  raw_data_matrix <- attr(x = seuratObj@assays[[assay]], which = slot)
  training_matrix <- as.data.frame(Matrix::t(as.matrix(raw_data_matrix)))
  training_matrix$celltype <- seuratObj@meta.data[,celltype_column]
  
  celltypes <- unique(training_matrix[,"celltype"])
  
  #Create output directories
  #Trained binary classifiers will be saved to /models
  #Parseable metrics .rds files will be saved to /metrics
  #Training data .rds files will be saved to /training_data
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    dir.create(paste0(output_dir,"/models"))
    dir.create(paste0(output_dir,"/metrics"))
    dir.create(paste0(output_dir, "/training_data"))
    print(paste0("Output directory: ", output_dir))
  } else {
    print(paste0("Output directory: ", output_dir))
  }
  

  
  #Iterate over celltypes and train a binary classifier for each celltype present in the celltype_column
  for (celltype in celltypes){
    temp_training_matrix <- training_matrix
    temp_training_matrix[,celltype] <- ifelse(training_matrix[,"celltype"]==celltype,1,0)
    names(temp_training_matrix) <- make.names(names(temp_training_matrix),unique = T)
    temp_model <- TrainModel(temp_training_matrix, celltype, hyperparameter_tuning = hyperparameter_tuning, learner = learner, inner_resampling = inner_resampling, outer_resampling = outer_resampling, inner_folds = inner_folds,inner_ratio = inner_ratio,  outer_folds = outer_folds, outer_ratio = outer_ratio, n_models = n_models, parallel = parallel)
    
    #trim the "celltype" column (leaving just the labeled varible celltype column as truth) and save the training matrix
    temp_training_matrix <- subset(temp_training_matrix, select = -celltype)
    saveRDS(temp_training_matrix, file = paste0(output_dir, "/training_data/", celltype, "_Training_Matrix.rds"))
   
    if(hyperparameter_tuning){
      #Save the trained model to the output directory
      saveRDS(temp_model$model, file = paste0(output_dir, "/models/", celltype, "_BinaryClassifier.rds"))
      saveRDS(temp_model$metrics, file = paste0(output_dir, "/metrics/", celltype, "_BinaryClassifier_Resampled.rds"))
    } else if (!hyperparameter_tuning){
      #Save the trained model to the output directory
      saveRDS(temp_model$model, file = paste0(output_dir, "/models/", celltype, "_BinaryClassifier.rds"))
      saveRDS(temp_model$metrics, file = paste0(output_dir, "/metrics/", celltype, "_BinaryClassifier_ConfusionMatrix.rds"))
    }
  }
  print("All models trained!")
  
  #Parse and print accuracy metrics from metrics rds files
  if(verbose){
    #Grab model names from model directory
    metrics_files <- list.files(paste0(output_dir, "/metrics"))
    for (metrics_file in metrics_files){
      print(metrics_file)
      .ParseMetricsFile(paste0(output_dir,"/metrics/",metrics_file))
    }
  }
}

#' @title Helper function for printing information about the metrics files
#'
#' @description Parses the .rds file saved during TrainAllModels and prints accuracy information
#' @param metrics_file Path to a metrics file written by TrainAllModels()
#' @export

.ParseMetricsFile<- function(metrics_file){
  if (grepl("ConfusionMatrix", metrics_file)){
    confusion <- readRDS(metrics_file)
    print(confusion$table)
    print(confusion$overall["Accuracy"])
  } else if (grepl("Resampled", metrics_file)){
    full_resampling <- readRDS(metrics_file)
    print(extract_inner_tuning_results(full_resampling))
    print(full_resampling$score())
    print(full_resampling$aggregate())
  } else{
    print("Unexpected metrics file filename. Please ensure the metrics files are generated by TrainAllModels or adhere to the naming convention if manually generated.")
    stop()
  }
}

#' @title Applies trained models to get celltype probabilities
#'
#' @description Applies trained models to get celltype probabilities.
#' @param seuratObj The Seurat Object to be updated
#' @param models_dir the directory containing the models produced by TrainAllModels
#' @param iterative Logical determining whether or not the data should be chunked to save memory or not. 
#' @param iterations Integer determing how many chunks the data should be split into. 
#' @export

PredictCellTypeProbability <- function(seuratObj, models_dir = "./classifiers/models/", iterative = F, iterations = 100){
  #Grab model names from model directory
  model_names <- list.files(models_dir)
  for (model_name in model_names){
    #Grab celltype name from the name of the model
    celltype <- strsplit(model_name,"_")[[1]][[1]]
    #load trained model
    print(paste("Classifying using",model_name))
    classifier <- readRDS(file = paste0(models_dir,model_name))
    #De-sparse and transpose seuratObj normalized data & make names unique
    gene_expression_matrix <- as.data.frame(Matrix::t(as.matrix(seuratObj@assays$RNA@data))) 
    names(gene_expression_matrix) <- make.names(names(gene_expression_matrix),unique = T)
    
    if (!iterative){
      #predict probabilities for current celltype
      probability_vector <- predict(classifier, newdata= gene_expression_matrix, predict_type = "prob")
      #append probabilities to seurat metadata
      seuratObj@meta.data[[paste0(celltype,"_probability")]] <- probability_vector[,1]
    } else if (iterative){
      #find intervals/chunks that split the data in the specified number of iterations
      intervals <- split(1:ncol(geneExpressionMatrix), sort(1:ncol(gene_expression_matrix) %% iterations))
      #instantiate a list to hold all of the probability vector
      probability_vector_list <- vector(mode = "list", length = iterations)
      i<-1
      for (interval in intervals){
        #Define interval/chunk limits
        start <- head(interval,1)
        end <- tail(interval,1)
        
        #report current iteration
        print(paste("Iteration", i , "of", iterations))
        #add vector of probabilities to list
        probability_vector_list[[i]] <- classifier$predict_newdata(model, gene_expression_matrix[,start:end])$prob[,1]
        i<-i+1
      }
      #convert from list of probabilities to vector
      probability_vector <- unlist(probability_vector_list)
      #append probabilities to seurat metadata
      seuratObj@meta.data[[,paste0(celltype,"_probability")]] <- probability_vector
    }
    
  }
  print("Classification finished!")
  return(seuratObj)   
}


#' @title Assigns celltype label based on probabilities in the metadata
#'
#' @description Assigns celltype label based on probabilities in the metadata
#' @param seuratObj The Seurat Object to be updated
#' @param minimum_probability The minimum probability for a confident cell type assignment
#' @param minimum_delta The minimum difference in probabilities necessary to call one celltype over another.
#' @export

AssignCellType <- function(seuratObj, minimum_probability = 0.5, minimum_delta = 0.25){
  #This grabs each column in the metadata with the suffix "_probability"
  probabilities_matrix <- seuratObj@meta.data[,grep("_probability",names(seuratObj@meta.data))]
  seuratObj@meta.data[,"Classifier_Consensus_Celltype"] <- "Unassigned"
  #Iterate over the cells in the seurat object
  for (cell in 1:nrow(probabilities_matrix)){
    #Find the name of the column with maximum probability and grab the celltype and store it as "top_label"
    max_probability_column <- which.max(probabilities_matrix[cell,])
    max_probability <- max(probabilities_matrix[cell,])
    top_label <- strsplit(names(max_probability_column),"_")[[1]][[1]]
    
    #Grab the second highest probability to compare with the max probability
    #the frustrating negated grepl expressions in these lines bypass the inability to evaluate a variable name as a column identifier in r (e.g. top_label would be parsed as literally the column name "top_label")
    second_highest_probability_column <- which.max(probabilities_matrix[cell, !grepl(top_label, names(probabilities_matrix))])
    second_highest_probability <- max(probabilities_matrix[cell, !grepl(top_label, names(probabilities_matrix))])
 
    #Check if the cell's highest probability classification exceeds the minimum probabilty set for a confident call.
    #Additionally check if the highest probability and second highest probability are at least minimum_delta apart. If not, assign Unknown.
    seuratObj@meta.data[cell,"Classifier_Consensus_Celltype"] <- ifelse( ((max_probability >= minimum_probability) & ((max_probability - second_highest_probability) > minimum_delta)), yes =  top_label , no = "Unknown")
    
  }
  return(seuratObj)
}

#' @title Interprets the feature importance of each model 
#'
#' @description Interprets the feature importance of each model using DALEX and feature importance permuation 
#' @param output_dir The output directory that TrainAllModels saved training data and models into
#' @param plot_type Argument to pass to model_parts(). Ratio or difference is recommended for large feature sets where the base model's AUC loss will be outside the plotting range.
#' @export

InterpretModels <- function(output_dir= "./classifiers", plot_type = "ratio"){
  #Sanitize output_dir
  if (endsWith(output_dir, "/")){
    print("Please remove the trailing slash in output_dir")
    stop()
  }
  
  #Iterate through models in the models directory
  model_names <- list.files(paste0(output_dir, "/models"))
  for (model_name in model_names){
    #Read in model
    model <- readRDS(paste0(output_dir, "/models/",model_name))
    #Grab celltype from model filename
    celltype <- strsplit(model_name,"_")[[1]][[1]]
    #Get associated training data
    training_data <- readRDS(paste0(output_dir, "/training_data/",celltype,"_Training_Matrix.rds"))
    #trim truth column
    data  <- training_data[, 1:(ncol(training_data)-1)]
    y <- training_data[, ncol(training_data)]

    #create explainer
    explainer <-DALEXtra::explain_mlr3(model    = model,
                             data     = data,
                             y        = y,
                             label    = celltype,
                             colorize = FALSE)
    
    #variable feature importance
    parts <- DALEX::model_parts(explainer, type = plot_type)
    print(plot(parts, max_vars=12, show_boxplots = FALSE))
    
  }
}
  
