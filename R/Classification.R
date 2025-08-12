#' @include Utils.R
#' @include CellTypist.R

#' @import mlr3
#' @import mlr3learners
#' @import ranger
#' @import ggplot2
#'

utils::globalVariables(
  names = c('Freq', 'Passed', 'Reference', 'Second', 'Top', 'TopLabel', 'booster', 'type', 'value', 'splitrule'),
  package = 'RIRA',
  add = TRUE
)

#' @title Creates a binary classifier to classify cells within a Seurat object
#'
#' @description Creates a binary classifier to classify cells
#' @param training_matrix A counts or data slot provided by TrainModelsFromSeurat
#' @param celltype The celltype (provided by TrainModelsFromSeurat) used as classifier's positive prediction
#' @param hyperparameter_tuning logical that determines whether or not hyperparameter tuning should be performed.
#' @param learner The mlr3 learner that should be used. Currently fixed to "classif.ranger" if hyperparameter tuning is FALSE. Otherwise, "classif.xgboost" and "classif.ranger" are supported.
#' @param inner_resampling The resampling strategy that is used for hyperparameter optimization. Holdout ("hout" or "holdout") and cross validation ("cv" or "cross-validation") are supported.
#' @param outer_resampling The resampling strategy that is used to determine overfitting. Holdout ("hout" or "holdout") and cross validation ("cv" or "cross-validation") are supported.
#' @param inner_folds The number of folds to be used for inner_resampling if cross-valdiation is performed.
#' @param outer_folds The number of folds to be used for outer_resampling if cross-valdiation is performed.
#' @param inner_ratio The ratio of training to testing data to be used for inner_resampling if holdout resampling is performed.
#' @param outer_ratio The ratio of training to testing data to be used for inner_resampling if holdout resampling is performed.
#' @param n_models The number of models to be trained during hyperparameter tuning. The model with the highest accuracy will be selected and returned.
#' @param n_cores If non-null, this number of workers will be used with future::plan
#' @export

TrainModel <- function(training_matrix, celltype, hyperparameter_tuning = F, learner = "classif.ranger", inner_resampling = "cv", outer_resampling = "cv", inner_folds = 4, inner_ratio = 0.8,  outer_folds = 3, outer_ratio = 0.8, n_models = 20, n_cores = NULL){
  set.seed(GetSeed())

  if (!is.null(n_cores)){
    future::plan("multisession", workers = n_cores)
  }

  #Fix gene names to conform with Seurat Object processing
  colnames(training_matrix) <- gsub(names(training_matrix), pattern = "-", replacement = ".")

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

  task <- mlr3::TaskClassif$new(classification.data, id = "CellTypeBinaryClassifier", target = "celltype_binary")

  if (!hyperparameter_tuning){
    #Use a ranger random forest tree with default parameters and holdout (80% train, 20% test) resampling
    #classification.data has column number = number_of_genes + 1, so the mtry argument uses the full gene matrix
    learner <- mlr3::lrn("classif.ranger", importance = "permutation", num.trees=500, mtry=(ncol(classification.data)-1), predict_type = "prob")
    set_threads(learner)
    train_set <- sample(task$nrow, 0.8 * task$nrow)
    test_set <- setdiff(seq_len(task$nrow), train_set)
    model <- learner$train(task, row_ids = train_set)$predict(task, row_ids=test_set)
    confusion <- caret::confusionMatrix(factor(model$response), factor(model$truth))

    return(list(model = learner, metrics = confusion))

  } else {
    #Set model-independent values for the autotuner
    measure <- msr("classif.ce")
    terminator <- mlr3verse::trm("evals", n_evals = n_models)

    #Define a tuning space 25% as large as the number of models
    #In the case of sensitive hyperparameters, resolution = 5 allows for a low/medium-low/medium/medium-high/high type parameter space
    tuner <- mlr3verse::tnr("grid_search", resolution = 5)

    #Define resampling method used for hyperparameter tuning
    if (inner_resampling == "cv" || inner_resampling == "cross-validation"){
      inner_resample <- rsmp("cv", folds = inner_folds)
    } else if (inner_resampling == "hout" || inner_resampling == "holdout"){
      inner_resample <- rsmp("holdout", ratio = inner_ratio)
    } else {
      stop("Unknown inner_resampling method provided. Please select one of cross-validation, cv, hout, or holdout")
    }

    #Define resampling method used to determine overfitting
    if (outer_resampling == "cv" || outer_resampling == "cross-validation"){
      outer_resample <- rsmp("cv", folds = outer_folds)
    } else if (outer_resampling == "hout" || outer_resampling == "holdout"){
      outer_resample <- rsmp("holdout", ratio = outer_ratio)
    } else {
      stop("Unknown outer_resampling method provided. Please select one of cross-validation, cv, hout, or holdout")
    }

    #Set learner and define parameter space
    if (learner == "classif.ranger"){
      #Define learner
      learner <- mlr3::lrn("classif.ranger", importance = "permutation", predict_type = "prob")

      #Define Ranger Hyperparameter Space (RandomBotv2)
      tune_ps <- mlr3verse::ps(
        num.trees = mlr3verse::p_int(lower = 10, upper = 2000),
        sample.fraction = mlr3verse::p_dbl(lower = 0.1, upper = 1),
        respect.unordered.factors = mlr3verse::p_fct(levels = c("ignore", "order", "partition")),
        min.node.size = mlr3verse::p_int(lower = 1, upper = 100),
        splitrule = mlr3verse::p_fct(levels = c("gini", "extratrees")),
        num.random.splits = mlr3verse::p_int(lower = 1, upper = 100, depends = splitrule == "extratrees")
        )
    } else if (learner == "classif.xgboost"){
      #Update task
      task <- mlr3::TaskClassif$new(classification.data, id = "CellTypeBinaryClassifier", target = "celltype_binary")
      #Define learner
      learner <- mlr3::lrn("classif.xgboost", predict_type = "prob")
      #Define XGBoost model's Hyperparameter Space (RandomBotv2)
      tune_ps <- mlr3verse::ps(
        booster = mlr3verse::p_fct(levels = c("gblinear", "gbtree", "dart")),
        nrounds = mlr3verse::p_int(lower = 2, upper = 8, trafo = function(x) as.integer(round(exp(x)))),
        eta = mlr3verse::p_dbl(lower = -4, upper = 0, trafo = function(x) 10^x),
        gamma = mlr3verse::p_dbl(lower = -5, upper = 1, trafo = function(x) 10^x),
        lambda = mlr3verse::p_dbl(lower = -4, upper = 3, trafo = function(x) 10^x),
        alpha = mlr3verse::p_dbl(lower = -4, upper = 3, trafo = function(x) 10^x),
        subsample = mlr3verse::p_dbl(lower = 0.1, upper = 1),
        max_depth = mlr3verse::p_int(lower = 1, upper = 15),
        min_child_weight = mlr3verse::p_dbl(lower = -1, upper = 0, trafo = function(x) 10^x),
        colsample_bytree = mlr3verse::p_dbl(lower = 0.1, upper = 1),
        colsample_bylevel = mlr3verse::p_dbl(lower = 0.1, upper = 1),
        rate_drop = mlr3verse::p_int(lower = 0, upper = 1, depends = booster == 'dart'),
        skip_drop = mlr3verse::p_int(lower = 0, upper = 1, depends = booster == 'dart')
      )
    }
  }

  #Define the autotuner using the parameter spaces and conditions defined above
  at <- mlr3tuning::AutoTuner$new(
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
#' @param n_cores If non-null, this number of workers will be used with future::plan
#' @param gene_list If non-null, the input count matrix will be subset to these features
#' @param gene_exclusion_list If non-null, the input count matrix will be subset to drop these features
#' @param verbose Whether or not to print the metrics data for each model after training.
#' @param min_cells_per_class If provided, any classes (and corresponding cells) with fewer than this many cells will be dropped from the training data
#' @export
TrainModelsFromSeurat <- function(seuratObj, celltype_column, assay = "RNA", slot = "data", output_dir = "./classifiers", hyperparameter_tuning = F, learner = "classif.ranger", inner_resampling = "cv", outer_resampling = "cv", inner_folds = 4, inner_ratio = 0.8,  outer_folds = 3, outer_ratio = 0.8, n_models = 20, n_cores = NULL, gene_list = NULL, gene_exclusion_list = NULL, verbose = TRUE, min_cells_per_class = 20){
  if (methods::missingArg(celltype_column)) {
    stop('Must provide the celltype_column argument')
  }

  if (!celltype_column %in% names(seuratObj@meta.data)) {
    stop(paste0('The column: ', celltype_column, ' is not present in the seurat object'))
  }

  if (endsWith(output_dir, "/")){
    output_dir <- gsub(output_dir, pattern = "/$", replacement = "")
  }

  if (!is.null(min_cells_per_class) && min_cells_per_class > 0) {
    seuratObj <- .DropLowCountClasses(seuratObj, celltype_column, min_cells_per_class)
  }

  #Read the raw data from a seurat object and parse into an mlr3-compatible labeled matrix
  raw_data_matrix <- attr(x = seuratObj@assays[[assay]], which = slot)
  if (!all(is.null(gene_list))) {
    gene_list <- ExpandGeneList(gene_list)
    if (!all(gene_list %in% rownames(raw_data_matrix))) {
      missing <- gene_list[!gene_list %in% rownames(raw_data_matrix)]
      stop(paste0('All features in gene_list must be present in the Seurat object features. Missing: ', paste0(missing, collapse = ',')))
    }

    raw_data_matrix <- raw_data_matrix[gene_list,]
  }

  if (!all(is.null(gene_exclusion_list))) {
    gene_exclusion_list <- ExpandGeneList(gene_exclusion_list)
    raw_data_matrix <- raw_data_matrix[!rownames(raw_data_matrix) %in% gene_exclusion_list,]
  }

  training_matrix <- as.data.frame(Matrix::t(as.matrix(raw_data_matrix)))
  training_matrix$celltype <- seuratObj@meta.data[,celltype_column]

  celltypes <- unique(training_matrix[,"celltype"])

  #Create output directories
  #Trained binary classifiers will be saved to /models
  #Parseable metrics .rds files will be saved to /metrics
  #Training data .rds files will be saved to /training_data
  for (fn in c(output_dir, paste0(output_dir,"/models"), paste0(output_dir,"/metrics"), paste0(output_dir, "/training_data"))) {
    if (!dir.exists(fn)){
      dir.create(fn)
    }
  }

  print(paste0("Output directory: ", output_dir))

  # TODO: rather than use cell types directly as file names, we should use make.names() or something to ensure they are valid and sane (i.e. 'CD8+ T cells')
  # TODO: rather than saving one RDS per classifier, would it make more sense to save a list of cellType -> classifier? This bundles everything into one file on disk?

  #Iterate over celltypes and train a binary classifier for each celltype present in the celltype_column
  print(paste0("Total cell types: ", length(celltypes)))
  for (celltype in celltypes){
    print(paste0("Training: ", celltype))
    temp_training_matrix <- training_matrix
    temp_training_matrix[,celltype] <- ifelse(training_matrix[,"celltype"]==celltype,1,0)

    if (sum(duplicated(names(temp_training_matrix))) > 0) {
      stop(paste0('Found duplicate names in the input matrix: ', paste0(names(temp_training_matrix)[duplicated(names(temp_training_matrix))], collapse = ',')))
    }
    names(temp_training_matrix) <- make.names(names(temp_training_matrix), unique = T)
    temp_model <- TrainModel(temp_training_matrix, celltype, hyperparameter_tuning = hyperparameter_tuning, learner = learner, inner_resampling = inner_resampling, outer_resampling = outer_resampling, inner_folds = inner_folds,inner_ratio = inner_ratio,  outer_folds = outer_folds, outer_ratio = outer_ratio, n_models = n_models, n_cores = n_cores)

    #trim the "celltype" column (leaving just the labeled varible celltype column as truth) and save the training matrix
    temp_training_matrix <- subset(temp_training_matrix, select = -celltype)
    saveRDS(temp_training_matrix, file = paste0(output_dir, "/training_data/", make.names(celltype), "_Training_Matrix.rds"))

    if (hyperparameter_tuning){
      #Save the trained model to the output directory
      saveRDS(temp_model$model, file = paste0(output_dir, "/models/", make.names(celltype), "_BinaryClassifier.rds"))
      saveRDS(temp_model$metrics, file = paste0(output_dir, "/metrics/", make.names(celltype), "_BinaryClassifier_Resampled.rds"))
    } else if (!hyperparameter_tuning){
      #Save the trained model to the output directory
      saveRDS(temp_model$model, file = paste0(output_dir, "/models/", make.names(celltype), "_BinaryClassifier.rds"))
      saveRDS(temp_model$metrics, file = paste0(output_dir, "/metrics/", make.names(celltype), "_BinaryClassifier_ConfusionMatrix.rds"))
    }
  }
  print("All models trained!")

  #Parse and print accuracy metrics from metrics rds files
  if (verbose){
    #Grab model names from model directory
    data <- NULL
    metrics_files <- list.files(paste0(output_dir, "/metrics"))
    for (metrics_file in metrics_files){
      toAdd <- .ParseMetricsFile(paste0(output_dir,"/metrics/",metrics_file))

      if (all(is.null(data))) {
        data <- toAdd
      } else {
        data <- rbind(data, toAdd)
      }
    }

    print(ggplot(data, aes(x = type, y = value, fill = type)) +
      geom_bar(stat = 'identity', color = 'black') +
      labs(x = 'Type', y = 'Accuracy') +
      egg::theme_presentation(base_size = 12) +
      ggtitle("Accuracy") +
      ylim(0,1) +
      theme(
        legend.position = 'none'
      )
    )
  }
}

.ParseMetricsFile <- function(metrics_file){
  data <- NULL
  if (grepl("ConfusionMatrix", metrics_file)){
    confusion <- readRDS(metrics_file)

    label <- gsub(basename(metrics_file), pattern = '_BinaryClassifier.rds', replacement = '')
    dat <- as.data.frame(confusion$table)
    dat2 <- as.data.frame(prop.table(confusion$table))

    P1 <- ggplot2::ggplot(dat2, aes(x = Prediction, y = Reference, z = Freq, fill = Freq)) +
      geom_tile() +
      geom_text(data = dat, mapping = aes(x = Prediction, y = Reference, label = Freq), inherit.aes = FALSE, size = 14) +
      egg::theme_presentation(base_size = 12) +
      ggtitle(label) +
      scale_fill_gradient2() +
      theme(
        legend.position = 'none'
      )

    P2 <- ggplot2::ggplot(dat2, aes(x = Prediction, y = Reference, z = Freq, fill = Freq)) +
      geom_tile() +
      geom_text(data = dat2, mapping = aes(x = Prediction, y = Reference, label = Freq), inherit.aes = FALSE, size = 14) +
      egg::theme_presentation(base_size = 12) +
      ggtitle(label) +
      scale_fill_gradient2() +
      theme(
        legend.position = 'none'
      )

    print(P1 + P2)

    toAdd <- data.frame(type = label, metric = "Accuracy", value = confusion$overall["Accuracy"])
    if (all(is.null(data))) {
      data <- toAdd
    } else {
      data <- rbind(data, toAdd)
    }
  } else if (grepl("Resampled", metrics_file)){
    full_resampling <- readRDS(metrics_file)
    #print(mlr3tuning::extract_inner_tuning_results(full_resampling))
    print(full_resampling$score())
    print(full_resampling$aggregate())
  } else{
    stop("Unexpected metrics file filename. Please ensure the metrics files are generated by TrainModelsFromSeurat or adhere to the naming convention if manually generated.")
  }

  return(data)
}


#' @title Applies a trained binary classifier to get per-cell probabilities.
#'
#' @description Applies a trained model to get per-cell probabilities.
#' @param seuratObj The Seurat Object to be updated
#' @param model Either the full filepath to a model RDS file, or the name of a built-in model.
#' @param fieldToClass A list mapping the target field name in the seurat object to the classifier level. The latter is either numeric, or the string label. For example: list('CD4_T' = 1, 'CD8_T' = 2))
#' @param batchSize To conserve memory, data will be chunked into batches of at most this many cells
#' @param assayName The assay holding gene expression data
#' @export
ScoreCellsWithSavedModel <- function(seuratObj, model, fieldToClass, batchSize = 20000, assayName = 'RNA') {
  classifier <- .ResolveModel(modelFile = model)

  #De-sparse and transpose seuratObj normalized data & make names unique
  gene_expression_matrix <- Matrix::t(Seurat::GetAssayData(seuratObj, assay = assayName, slot = "data"))

  # NOTE: makeNames() will convert hyphen to period, and also prefix genes with numeric starts, like 7SK.2 -> X7SK.2
  colnames(gene_expression_matrix) <- make.names(colnames(gene_expression_matrix))

  # TODO: Is there is more universal way to get the model to report its features?
  # See: https://mlr3.mlr-org.com/reference/LearnerClassif.html
  modelFeats <- NULL
  if ('importance' %in% classifier$properties) {
    modelFeats <- names(classifier$importance())
  } else if (!is.null(classifier$model) && grepl(x = classifier$model$TypeDetail, pattern = 'logistic regression')) {
    modelFeats <- colnames(classifier$model$W)
    toRemove <- names(classifier$param_set$levels)
    modelFeats <- modelFeats[!modelFeats %in% toRemove]
  }

  if (!all(is.null(modelFeats))){
    missing <- modelFeats[!modelFeats %in% colnames(gene_expression_matrix)]
    if (length(missing) > 0) {
      stop(paste0('The following features are used in the model and missing from the input: ', paste0(sort(missing), collapse = ',')))
    }

    # Subset input data to match model:
    toDrop <- !(colnames(gene_expression_matrix) %in% modelFeats)
    if (sum(toDrop) > 0) {
      print(paste0('Dropping features not present in model: ', sum(toDrop), ' of ', ncol(gene_expression_matrix)))
      gene_expression_matrix <- gene_expression_matrix[,!toDrop, drop = FALSE]
    }
  } else {
    warning(paste0('Unable to infer features from model, type: ', classifier$model$TypeDetail))
  }

  print(paste0('Features shared between gene matrix and model: ', ncol(gene_expression_matrix)))

  nBatches <- ifelse(is.na(batchSize), yes = 1, no = ceiling(nrow(gene_expression_matrix) / batchSize))
  probability_vectors <- list()
  for (batchIdx in 1:nBatches){
    start <- 1 + ((batchIdx-1) * batchSize)
    end <- min((batchIdx * batchSize), nrow(gene_expression_matrix))
    print(paste0("Iteration ", batchIdx, " of ", nBatches, ", (", start, "-", end, ")"))

    #columns are named '0','1', so the first column is '0' and the second column is '1'
    dat <- stats::predict(classifier, newdata = data.frame(gene_expression_matrix[start:end,]), predict_type = 'prob')

    for (fieldName in names(fieldToClass)) {
      idx <- fieldToClass[[fieldName]]
      if (is.na(as.numeric(idx))) {
        # Try to resolve from model:
        if (idx %in% levels(classifier$model$ClassNames)) {
          idx <- which(levels(classifier$model$ClassNames) == idx)
        } else {
          stop(paste0('Unknown class: ', idx))
        }
      } else {
        idx <- as.numeric(idx)
      }

      if (batchIdx == 1) {
        probability_vectors[[fieldName]] <- dat[,idx]
      } else {
        probability_vectors[[fieldName]] <- c(probability_vectors[[fieldName]], dat[,idx])
      }
    }
  }

  for (fieldName in names(probability_vectors)) {
    probability_vector <- probability_vectors[[fieldName]]

    if (length(probability_vector) != ncol(seuratObj)) {
      stop(paste0('Error calculating probability_vector. Length was: ', length(probability_vector)))
    }

    #append probabilities to seurat metadata
    seuratObj@meta.data[[fieldName]] <- probability_vector

    if (length(names(seuratObj@reductions)) > 0) {
      print(Seurat::FeaturePlot(seuratObj, features = fieldName))
    }
  }

  return(seuratObj)
}


#' @title Applies trained models to get celltype probabilities
#'
#' @description Applies trained models to get celltype probabilities.
#' @param seuratObj The Seurat Object to be updated
#' @param models A named vector of models, where the values are the modelName (for built-in model), or filePath to an RDS file. The names of the vector should be the cell-type label for cells scored positive by the classifier.
#' @param fieldName The name of the metadata column to store the result
#' @param batchSize To conserve memory, data will be chunked into batches of at most this many cells
#' @param assayName The assay holding gene expression data
#' @param minimum_probability The minimum probability for a confident cell type assignment
#' @param minimum_delta The minimum difference in probabilities necessary to call one celltype over another.
#' @export
PredictCellTypeProbability <- function(seuratObj, models, fieldName = 'RIRA_Consensus', batchSize = 20000, assayName= "RNA", minimum_probability = 0.5, minimum_delta = 0.25){
  fieldNames <- c()
  for (modelName in names(models)) {
    print(paste0('Scoring with model: ', modelName))
    probColName <- paste0(modelName, '_probability')
    fieldNames <- c(fieldNames, probColName)
    fieldToClass <- list()
    fieldToClass[[probColName]] <- 2 # Assume binary classifier for now
    seuratObj <- ScoreCellsWithSavedModel(seuratObj, model = models[[modelName]], fieldToClass = fieldToClass, batchSize = batchSize, assayName = assayName)
  }

  seuratObj <- AssignCellType(seuratObj, probabilityColumns = fieldNames, fieldName = fieldName, minimum_probability = minimum_probability, minimum_delta = minimum_delta)

  return(seuratObj)
}


#' @title Assigns celltype label based on probabilities in the metadata
#'
#' @description Assigns celltype label based on probabilities in the metadata
#' @param seuratObj The Seurat Object to be updated
#' @param probabilityColumns The set of columns containing probabilities for each classifier to include
#' @param fieldName The name of the metadata column to store the result
#' @param minimum_probability The minimum probability for a confident cell type assignment
#' @param minimum_delta The minimum difference in probabilities necessary to call one celltype over another.
#' @export
AssignCellType <- function(seuratObj, probabilityColumns, fieldName = 'RIRA_Consensus', minimum_probability = 0.5, minimum_delta = 0.25){
  probabilities_matrix <- seuratObj@meta.data[,probabilityColumns, drop = F]
  if (ncol(probabilities_matrix) == 0) {
    stop('Unable to find cell type probability columns!')
  }

  seuratObj@meta.data[,fieldName] <- "Unassigned"
  #Iterate over the cells in the seurat object

  toPlot <- NULL
  for (cell in 1:nrow(probabilities_matrix)){
    #Find the name of the column with maximum probability and grab the celltype and store it as "top_label"
    max_probability_column <- which.max(probabilities_matrix[cell,])
    max_probability <- max(probabilities_matrix[cell,])
    top_label <- strsplit(names(max_probability_column),"_")[[1]][[1]]
    second_highest_probability <- max(probabilities_matrix[cell, names(probabilities_matrix) != names(max_probability_column)])

    #Check if the cell's highest probability classification exceeds the minimum probabilty set for a confident call.
    #Additionally check if the highest probability and second highest probability are at least minimum_delta apart. If not, assign Unknown.
    seuratObj@meta.data[cell,fieldName] <- ifelse( ((max_probability >= minimum_probability) & ((max_probability - second_highest_probability) > minimum_delta)), yes =  top_label , no = "Unknown")

    passed <- seuratObj@meta.data[cell,fieldName] == top_label
    toAdd <- data.frame(CellBarcode = colnames(seuratObj)[cell], Top = max_probability, TopLabel = top_label, Second = second_highest_probability, Passed = passed)
    if (all(is.null(toPlot))) {
      toPlot <- toAdd
    } else {
      toPlot <- rbind(toPlot, toAdd)
    }
  }

  minVal <- min(c(toPlot$Top, toPlot$Second)) * 0.9  # provide consistent x/y limits
  P1 <- ggplot2::ggplot(toPlot, aes(x = Top, y = Second, color = TopLabel, shape = Passed)) +
    geom_point() +
    ggtitle("Cell Type Probabilities") +
    egg::theme_presentation(base_size = 10) +
    ylim(minVal, 1) +
    xlim(minVal, 1) +
    labs(x = 'Highest Probability', y = 'Second Highest Probability', color = 'Top Label', shape = 'Passed?')

  print(P1)

  return(seuratObj)
}

.ResolveModel <- function(modelFile) {
  if (file.exists(modelFile)) {
    return(readRDS(file = modelFile))
  }

  savedModel <- system.file(paste0("models/", modelFile, ".rds"), package = "RIRA")
  if (!file.exists(savedModel)) {
    stop(paste0('Unable to find model: ', modelFile))
  }

  return(readRDS(file = savedModel))
}

#' @title Interprets the feature importance of each model
#'
#' @description Interprets the feature importance of each model using DALEX and feature importance permuation
#' @param output_dir The output directory that TrainAllModels saved training data and models into
#' @param plot_type Argument to pass to model_parts(). Ratio or difference is recommended for large feature sets where the base model's AUC loss will be outside the plotting range.
#' @export
InterpretModels <- function(output_dir= "./classifiers", plot_type = "ratio"){
  if (endsWith(output_dir, "/")){
    output_dir <- gsub(output_dir, pattern = "/$", replacement = "")
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
    explainer <- DALEXtra::explain_mlr3(model    = model,
                                        data     = data,
                                        y        = y,
                                        label    = celltype,
                                        colorize = FALSE)

    #variable feature importance
    parts <- DALEX::model_parts(explainer, type = plot_type)
    print(plot(parts, max_vars=12, show_boxplots = FALSE))

  }
}
#' @title Predicts T cell activation using sPLS derived components and a trained logistic model on transformed variates
#' @description Predicts T cell activation using a trained model
#' @param seuratObj The Seurat Object to be updated
#' @param model The trained sPLS model to use for prediction. This can be a file path to an RDS file, or a built-in model name.
#' @return A Seurat object with the sPLS scores and predicted probabilities added to the metadata
#' @export

PredictTcellActivation <- function(seuratObj, model = NULL) {
  #check & sanitize model
  if (is.null(model)) {
    #default method
    print("No model provided, using RIRA's built-in T Cell Activation model.")
    modelFile <- system.file("models/ActivatedTCell4ClassModel_v1.rds", package = "RIRA")
    model <- readRDS(modelFile)
  } else if (grepl(pattern = "\\.rds$", x = model, fixed = TRUE) && file.exists(model)) {
    #user provides .rds
    print(paste0("Using model from file: ", model))
    model <- readRDS(model)
  } else if (!.CanPredict(model)) {
    #user provides some other kind of model object, but it can't predict using stats::predict
    stop(paste0("Provided model: ", model, " does not have a detectable predict method. Please provide a valid model or file path to an RDS file containing a trained model."))
  } else {
    stop("Model must be one of: NULL (built-in), a file path to an RDS file, built-in model name. If you want to use a custom model, please ensure it is compatible with stats::predict().\nIf you want to use the built-in model, please provide NULL as the model argument.\n")
  }

  modelCoefs <- colnames(stats::coef(model))
  if (!all(paste0("comp",1:6) %in% modelCoefs)) {
    if ('nnet' %in% class(model)) {
      modelCoefs <- colnames(nnet:::coef.multinom(model))
    }

    if (!all(paste0("comp",1:6) %in% modelCoefs)) {
      stop("Model does not contain the expected components. Please ensure the model's features are conformant with this prediction method.\nExpected components: comp1, comp2, comp3, comp4, comp5, comp6. Found: " + paste0(colnames(stats::coef(model)), collapse = ','))
    }
  }

  #define components & score
  comps <- paste0("PLS_Score_", seq_len(6))
  for (i in comps){
    seuratObj <- ScoreUsingSavedComponent(seuratObj, componentOrName = i, fieldName = i, layer = "scale.data")
  }
  #construct prediction for the model
  newdata <- Seurat::FetchData(seuratObj, vars = paste0("PLS_Score_", seq_len(6)))
  colnames(newdata) <- paste0("comp", seq_len(6))

  if (!all(rownames(newdata) == colnames(seuratObj))) {
    stop("Internal Error: Cell names in FetchData() do not match Seurat cell names.")
  }
  #predictions/scoring
  prob_df <- stats::predict(model, newdata, type = "prob")
  class_vec <- stats::predict(model, newdata, type = "class")

  colnames(prob_df) <- paste0("sPLS_prob_", colnames(prob_df))
  class_df <- data.frame(
    sPLS_class = as.character(class_vec),
    row.names = rownames(prob_df),
    stringsAsFactors = FALSE
  )
  #add back to seurat & return
  seuratObj <- Seurat::AddMetaData(seuratObj, metadata = prob_df)
  seuratObj <- Seurat::AddMetaData(seuratObj, metadata = class_df)
  return(seuratObj)
}

#basic predict wrapper to check if the model can be used with stats::predict()
.CanPredict <- function(model, newdata = NULL) {
  tryCatch({
    if (is.null(newdata)) {
      #naive predict call, no newdata provided
      result <- predict(model)
    } else {
      result <- predict(model, newdata = newdata)
    }
    return(TRUE)
  }, error = function(e) {
    error_msg <- tolower(conditionMessage(e))
    #catch common errors that suggest the model is compatible with stats::predict()
    if (grepl("newdata|argument.*missing|unused argument", error_msg)) {
      return(TRUE)
    }
    #catch errors that suggest the predict method is not defined for the model
    if (grepl("no applicable method|could not find function", error_msg)) {
      return(FALSE)
    }
    #for other errors, assume method exists but there's a usage issue
    return(TRUE)
  })
}
