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
#' @param model The trained sPLSDA model to use for prediction. This can be a file path to an RDS file, or a built-in model name.
#' @param modelList A list of trained sPLSDA models to use for prediction. This can be a list of file paths to RDS files, or built-in model names.
#' @import nnet
#' @return A Seurat object with the sPLSDA scores and predicted probabilities added to the metadata
#' @export

PredictTcellActivation <- function(seuratObj, model = NULL, modelList = NULL) {
  #################
  ### Sanitize  ###
  #################
  #check & sanitize model
  #track model paths and versions
  modelPaths <- list()
  modelVersions <- list()
  defaults <- FALSE
  
  if (is.null(model) && is.null(modelList)) {
    #default method
    print("No model provided, using RIRA's built-in CD8 and CD4 T Cell Activation models.")
    modelFiles <- list(
      CD8 = "models/ActivatedCD8TCell4ClassModel_v1.rds",
      CD4 = "models/ActivatedCD4TCell4ClassModel_v1.rds", 
      General = "models/GeneralActivatedTCellModel_v3.rds"
    )
    
    modelList <- list()
    for (modelName in names(modelFiles)) {
      modelPath <- system.file(modelFiles[[modelName]], package = "RIRA")
      modelList[[modelName]] <- readRDS(modelPath)
      modelPaths[[modelName]] <- modelFiles[[modelName]]
      #extract version from filename (e.g., _v1.rds -> v1)
      version <- gsub(".*_(v\\d+)\\.rds$", "\\1", modelFiles[[modelName]])
      modelVersions[[modelName]] <- version
    }
    #flag that the defaults were used
    defaults <- TRUE
  } else if (!is.null(model) && grepl(pattern = "\\.rds$", x = model, fixed = TRUE) && file.exists(model)) {
    #user provides .rds
    print(paste0("Using model from file: ", model))
    modelList <- list(Custom = readRDS(model))
    modelPaths[["Custom"]] <- model
    #try to extract version from filename, default to "custom" if no version found
    version <- gsub(".*_(v\\d+)\\.rds$", "\\1", basename(model))
    modelVersions[["Custom"]] <- ifelse(grepl("^v\\d+$", version), version, "custom")
  } else if (!is.null(model)) {
    stop("Model must be one of: NULL (built-ins), a file path to an RDS file. If you want to use a custom model, please ensure it is compatible with stats::predict().\nIf you want to use the built-in models, please provide NULL as the model argument.\n")
  }
  #check that models can predict
  for (modelName in names(modelList)) {
    modelObj <- modelList[[modelName]]
    if (!.CanPredict(modelObj)) {
      #user provides some other kind of model object, but it can't predict using stats::predict
      if (length(modelList) == 1) { 
        stop(paste0("Provided model does not have a detectable predict method. Please provide a valid model or file path to an RDS file containing a trained model."))
      } else {
        stop(paste0("Model '", modelName, "' does not have a detectable predict method. Please provide valid models or file paths to RDS files containing trained models."))
      }
    }
  }
  
  #######################################
  ### Default Model Component Scoring ###
  #######################################
  #if the defaults were used, we ensure the component scores are present and valid
  if (defaults) {
    #check if component scores already exist in metadata and score if missing or contain NAs
    metadataNames <- colnames(seuratObj@meta.data)
    for (modelName in names(modelList)) {
      #determine number of components for this model
      nComponents <- .GetModelComponentCount(modelName)
      modelVersion <- modelVersions[[modelName]]
      print(paste0("Model '", modelName, "' (", modelVersion, ") uses ", nComponents, " components"))
      
      needsScoring <- FALSE
      #check if all component scores exist for this model
      for (i in seq_len(nComponents)){
        fieldName <- paste0(modelName, "_Activation_sPLSDA_Score_", i, "_", modelVersion)
        if (!fieldName %in% metadataNames) {
          needsScoring <- TRUE
          print(paste0("Component score '", fieldName, "' not found in metadata. Will calculate."))
          break
        } else {
          #check for NAs in existing column
          if (any(is.na(seuratObj@meta.data[[fieldName]]))) {
            needsScoring <- TRUE
            print(paste0("Component score '", fieldName, "' contains NA values. Will recalculate."))
            break
          }
        }
      }
      #score components if needed
      if (needsScoring) {
        print(paste0("Scoring sPLSDA components for ", modelName, " model (", modelVersion, ")..."))
        componentPrefix <- paste0(modelName, "_Activation_sPLSDA_component")
        for (i in seq_len(nComponents)){
          componentName <- paste0(componentPrefix, i)
          fieldName <- paste0(modelName, "_Activation_sPLSDA_Score_", i, "_", modelVersion)
          seuratObj <- ScoreUsingSavedComponent(seuratObj, componentOrName = componentName, fieldName = fieldName, layer = "scale.data")
        }
      }
    }
  }
  
  ###############
  ### Scoring ###
  ###############
  #iterate models
  for (modelIdx in seq_along(modelList)) {
    modelName <- names(modelList)[modelIdx]
    modelObj <- modelList[[modelIdx]]
    modelVersion <- modelVersions[[modelName]]
    
    print(paste0("Processing model: ", modelName, " (", modelVersion, ")"))
    
    #determine expected number of components from model coefficients
    modelCoefs <- colnames(stats::coef(modelObj))
    if (!'nnet' %in% class(modelObj) || is.null(modelCoefs) || !any(grepl("^comp", modelCoefs))) {
      if ('nnet' %in% class(modelObj)) {
        modelCoefs <- colnames(stats::coef(modelObj))
      }
    }
    
    #extract component numbers from coefficient names
    compNames <- modelCoefs[grepl("^comp[0-9]+$", modelCoefs)]
    if (length(compNames) == 0) {
      stop("Model does not contain component features (comp1, comp2, etc.). Please ensure the model's features are conformant with this prediction method.\nFound coefficients: ", paste0(dplyr::coalesce(modelCoefs, 'NULL'), collapse = ','))
    }
    
    #determine number of components from the highest numbered component
    compNumbers <- as.numeric(gsub("comp", "", compNames))
    nComponents <- max(compNumbers)
    expectedComps <- paste0("comp", seq_len(nComponents))
    
    if (!all(expectedComps %in% modelCoefs)) {
      stop("Model does not contain all expected components. Please ensure the model's features are conformant with this prediction method.\nExpected components: ", paste0(expectedComps, collapse = ', '), ". Found: ", paste0(compNames, collapse = ','))
    }
    
    #construct prediction for the model, forcing conformance in component names
    #look for component scores with version number
    scoreVars <- paste0(modelName, "_Activation_sPLSDA_Score_", seq_len(nComponents), "_", modelVersion)
    newdata <- Seurat::FetchData(seuratObj, vars = scoreVars)
    colnames(newdata) <- paste0("comp", seq_len(nComponents))
    
    if (!all(rownames(newdata) == colnames(seuratObj))) {
      stop("Internal Error: Cell names in FetchData() do not match Seurat cell names.")
    }
    
    #predictions/scoring
    prob_df <- stats::predict(modelObj, newdata, type = "prob")
    class_vec <- stats::predict(modelObj, newdata, type = "class")
    
    #include version in column names
    colnames(prob_df) <- paste0(modelName, "_sPLS_prob_", colnames(prob_df), "_", modelVersion)
    class_df <- data.frame(
      sPLS_class = as.character(class_vec),
      row.names = rownames(prob_df),
      stringsAsFactors = FALSE
    )
    colnames(class_df) <- paste0(modelName, "_sPLS_class_", modelVersion)
    
    #add predictions back to Seurat object
    seuratObj <- Seurat::AddMetaData(seuratObj, metadata = prob_df)
    seuratObj <- Seurat::AddMetaData(seuratObj, metadata = class_df)
    
    if (length(names(seuratObj@reductions)) > 0) {
      print(Seurat::DimPlot(seuratObj, group.by = paste0(modelName, "_sPLS_class_", modelVersion)))
    }
  }

  seuratObj <- CombineTcellActivationClasses(seuratObj, classMapping = GetActivationClassMapping('TcellActivation.Basic'), outputFieldName = 'sPLS_TCR_General_v3')
  if (any(is.na(seuratObj$sPLS_TCR_General_v3))) {
    stop('There were NA values for seuratObj$sPLS_TCR_General_v3')
  }

  seuratObj$Is_TCR_Stimulated <- seuratObj$sPLS_TCR_General_v3 %in% c('General_Combined_prob_Th1-Tc1-like_v3', 'General_Combined_prob_Th2-Th17-like_v3')
  if (length(names(seuratObj@reductions)) > 0) {
    print(DimPlot(seuratObj, group.by = 'Is_TCR_Stimulated'))
  }

  return(seuratObj)
}

#helper function to determine the number of components for a built-in model
.GetModelComponentCount <- function(modelName) {
  #list all component files for this model
  componentPattern <- paste0(modelName, "_Activation_sPLSDA_component")
  componentFiles <- list.files(system.file("components", package = "RIRA"), 
                               pattern = paste0("^", componentPattern, "[0-9]+\\.tsv$"),
                               full.names = FALSE)
  
  if (length(componentFiles) == 0) {
    stop(paste0("No component files found for model '", modelName, "'. Expected files matching pattern: ", componentPattern, "N.tsv"))
  }
  
  #extract component numbers from filenames
  compNumbers <- as.numeric(gsub(paste0(componentPattern, "(\\d+)\\.tsv"), "\\1", componentFiles))
  
  #return the highest component number
  return(max(compNumbers))
}

#basic predict wrapper to check if the model can be used with stats::predict()
.CanPredict <- function(model, newdata = NULL) {
  tryCatch({
    if (is.null(newdata)) {
      #naive predict call, no newdata provided
      result <- stats::predict(model)
    } else {
      result <- stats::predict(model, newdata = newdata)
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

#' @title Combine T cell activation classes using custom logic
#' @description Takes a Seurat object with T cell activation predictions (from PredictTcellActivation) 
#' and combines classes based on user-defined logic. This is useful for collapsing fine-grained 
#' activation states into broader categories.
#' @param seuratObj The Seurat Object containing T cell activation predictions
#' @param modelName The name of the model whose predictions should be combined. Default is "General".
#' @param modelVersion The version of the model. Default is "v3" for the General model.
#' @param classMapping A named list where names are the new combined class labels and values are 
#' character vectors of the original classes to combine. You can also fetch predefined mappings
#' using GetActivationClassMapping(name). For example: GetActivationClassMapping('TcellActivation.Basic')
#' or a manual list such as list("Activated" = c("Early_Activated", "Late_Activated"), "Resting" = c("Naive", "Memory")).
#' @param outputFieldName The name of the metadata column to store the combined classifications. 
#' If NULL, will use "<modelName>_Combined_Class_<modelVersion>".
#' @param probabilityAggregation How to aggregate probabilities for combined classes. Options are:
#' "sum" (default) - sum probabilities of constituent classes,
#' "max" - take maximum probability among constituent classes,
#' "mean" - take mean of probabilities of constituent classes.
#' @param relabelOrRecall Strategy to set the combined class column. Options:
#' "relabel" uses the original class-to-mapping lookup, labeling unmapped as "Unmapped".
#' "recall" (default) uses the newly aggregated combined probabilities to assign classes.
#' @param recallMethod When relabelOrRecall = "recall", how to recall combined classes.
#' Options: "max" (default) choose the highest combined probability, or "threshold" which
#' assigns the first combined class with probability >= recallProbabilityThreshold, otherwise "Uncalled".
#' @param recallProbabilityThreshold Threshold used when recallMethod = "threshold". Default: 0.5.
#' @return A Seurat object with the combined class assignments added to metadata
#' @seealso GetActivationClassMapping, PredictTcellActivation
#' @export
#' @examples
#' \dontrun{
#' # Default usage with all parameters explicit
#' mapping <- GetActivationClassMapping('TcellActivation.Basic')
#' seuratObj <- CombineTcellActivationClasses(
#'   seuratObj,
#'   modelName = "General",
#'   modelVersion = "v3",
#'   classMapping = mapping,
#'   outputFieldName = NULL,
#'   probabilityAggregation = "sum",
#'   relabelOrRecall = "recall",
#'   recallMethod = "max",
#'   recallProbabilityThreshold = 0.5
#' )
#'
#' # Using a pre-registered mapping
#' mapping <- GetActivationClassMapping('TcellActivation.Basic')
#' seuratObj <- CombineTcellActivationClasses(
#'   seuratObj,
#'   classMapping = mapping
#' )
#'
#' # Combine activation states into broader categories
#' seuratObj <- CombineTcellActivationClasses(
#'   seuratObj,
#'   classMapping = list(
#'     "Th1" = c("Th1_MIP1B.neg_CD137.neg", "Th1_MIP1B.neg_CD137.pos", "Th1_MIP1B.pos"),
#'     "Th17" = c("Th17"), 
#'     "Bystander Activated" = c("NonSpecificActivated_L1", "NonSpecificActivated_L2"), 
#'     "Bulk Tissue T cell" = c("Uncultured"), 
#'     "Naive T cell" = c("Cultured_Bystander_NoBFA", "Cultured_Bystander_BFA")
#'   )
#' )
#' 
#' # Using Guidelines for T cell nomenclature (https://www.nature.com/articles/s41577-025-01238-2)
#' # with a slight augmentation that antigen specific T cells are "star" for antigens, rather than plus or zero.
#' seuratObj <- CombineTcellActivationClasses(
#'   seuratObj,
#'   classMapping = list(
#'     "CD4posOrCD8pos_Th1_U_B_A_t_star" = c("Th1_MIP1B.neg_CD137.neg", "Th1_MIP1B.neg_CD137.pos", "Th1_MIP1B.pos"),
#'     "CD4pos_Th17_U_B_A_t_star" = c("Th17"), 
#'     "CD4posOrCD8pos_Th1_U_B_A_t_O" = c("NonSpecificActivated_L1", "NonSpecificActivated_L2"), 
#'     "CD4posOrCD8pos_T_U_R_N_t_O" = c("Uncultured"), 
#'     "CD4posOrCD8pos_T_U_B_N_t_O" = c("Cultured_Bystander_NoBFA", "Cultured_Bystander_BFA")
#'   )
#' )
#' }

CombineTcellActivationClasses <- function(seuratObj, 
                                          modelName = "General", 
                                          modelVersion = "v3",
                                          classMapping,
                                          outputFieldName = NULL,
                                          probabilityAggregation = "sum",
                                          relabelOrRecall = "recall",
                                          recallMethod = "max",
                                          recallProbabilityThreshold = 0.5) {
  #################
  ### Sanitize  ###
  #################
  if (missing(classMapping) || is.null(classMapping)) {
    stop("classMapping must be provided. It should be a named list mapping new class names to vectors of original class names.")
  }
  
  if (!is.list(classMapping) || is.null(names(classMapping))) {
    stop("classMapping must be a named list.")
  }
  
  #validate that all names in classMapping are non-empty strings and not the literal 'NULL'
  nameVec <- names(classMapping)
  if (any(nameVec == "" | is.na(nameVec))) {
    stop("All elements in classMapping must have non-empty names.")
  }
  if (any(tolower(nameVec) == "null")) {
    stop("classMapping names cannot be 'NULL'. Please provide meaningful combined class names.")
  }
  
  #validate that all values in classMapping are character vectors
  for (newClassName in names(classMapping)) {
    if (!is.character(classMapping[[newClassName]])) {
      stop(paste0("Value for '", newClassName, "' must be a character vector of original class names."))
    }
    if (length(classMapping[[newClassName]]) == 0) {
      stop(paste0("Value for '", newClassName, "' cannot be an empty vector. Please provide at least one original class name."))
    }
    if (any(is.na(classMapping[[newClassName]]) | classMapping[[newClassName]] == "")) {
      stop(paste0("Value for '", newClassName, "' contains empty or NA class names. All class names must be valid strings."))
    }
  }
  
  #check for duplicate classes across different mappings
  allMappedClasses <- unlist(classMapping, use.names = FALSE)
  duplicateClasses <- allMappedClasses[duplicated(allMappedClasses)]
  if (length(duplicateClasses) > 0) {
    stop(paste0("The following classes appear in multiple mappings: ", 
                paste0(unique(duplicateClasses), collapse = ", "), 
                ". Each original class can only be mapped to one combined class."))
  }
  
  if (!probabilityAggregation %in% c("sum", "max", "mean")) {
    stop("probabilityAggregation must be one of: 'sum', 'max', 'mean'")
  }

  if (!relabelOrRecall %in% c("relabel", "recall")) {
    stop("relabelOrRecall must be one of: 'relabel', 'recall'")
  }
  if (!recallMethod %in% c("max", "threshold")) {
    stop("recallMethod must be one of: 'max', 'threshold'")
  }
  if (!is.numeric(recallProbabilityThreshold) || length(recallProbabilityThreshold) != 1 || is.na(recallProbabilityThreshold)) {
    stop("recallProbabilityThreshold must be a single numeric value")
  }
  if (recallProbabilityThreshold < 0 || recallProbabilityThreshold > 1) {
    stop("recallProbabilityThreshold must be in [0, 1]")
  }
  
  #check if the model predictions exist in the metadata
  classFieldName <- paste0(modelName, "_sPLS_class_", modelVersion)
  if (!classFieldName %in% colnames(seuratObj@meta.data)) {
    stop(paste0("Class predictions '", classFieldName, "' not found in metadata. Please run PredictTcellActivation() first."))
  }
  
  #get all original classes from the model
  originalClasses <- unique(seuratObj@meta.data[[classFieldName]])
  
  #check that all classes in classMapping exist in the original predictions
  mappedClasses <- allMappedClasses
  missingClasses <- setdiff(mappedClasses, originalClasses)
  if (length(missingClasses) > 0) {
    warning(paste0("The following classes in classMapping were not found in the original predictions and will be ignored for probability aggregation: ", 
                   paste0(missingClasses, collapse = ", ")))
  }
  
  #warn if some original classes are not mapped
  unmappedClasses <- setdiff(originalClasses, mappedClasses)
  if (length(unmappedClasses) > 0) {
    warning(paste0("The following classes were present in predictions but not accounted for in classMapping and will be labeled as 'Unmapped': ",
                   paste0(unmappedClasses, collapse = ", ")))
  }
  
  #set output field name
  if (is.null(outputFieldName)) {
    outputFieldName <- paste0(modelName, "_Combined_Class_", modelVersion)
  }
  
  ########################
  ### Combine Classes  ###
  ########################
  #initialize combined class vector
  combinedClasses <- character(ncol(seuratObj))
  names(combinedClasses) <- colnames(seuratObj)
  
  #initialize combined probability matrix
  combinedProbs <- matrix(0, nrow = ncol(seuratObj), ncol = length(classMapping))
  colnames(combinedProbs) <- names(classMapping)
  rownames(combinedProbs) <- colnames(seuratObj)
  
  #get all probability columns for this model
  probPattern <- paste0("^", modelName, "_sPLS_prob_.*_", modelVersion, "$")
  probCols <- grep(probPattern, colnames(seuratObj@meta.data), value = TRUE)
  
  if (length(probCols) == 0) {
    warning(paste0("No probability columns found for model '", modelName, "' version '", modelVersion, 
                   "'. Combined probabilities will not be calculated."))
  }
  
  # if relabelOrRecall == "relabel", map each cell to its combined class using original predictions
  if (relabelOrRecall == "relabel") {
    for (cellIdx in seq_len(ncol(seuratObj))) {
      originalClass <- seuratObj@meta.data[cellIdx, classFieldName]
      #find which combined class this original class belongs to
      combinedClass <- "Unmapped"
      for (newClassName in names(classMapping)) {
        if (originalClass %in% classMapping[[newClassName]]) {
          combinedClass <- newClassName
          break
        }
      }
      combinedClasses[cellIdx] <- combinedClass
    }
  }
  
  #aggregate probabilities for combined classes
  if (length(probCols) > 0) {
    for (newClassName in names(classMapping)) {
      constituentClasses <- classMapping[[newClassName]]
      
      #find probability columns for constituent classes
      constituentProbCols <- c()
      for (constituent in constituentClasses) {
        probCol <- paste0(modelName, "_sPLS_prob_", constituent, "_", modelVersion)
        if (probCol %in% probCols) {
          constituentProbCols <- c(constituentProbCols, probCol)
        }
      }
      
      if (length(constituentProbCols) > 0) {
        #aggregate probabilities based on specified method
        if (probabilityAggregation == "sum") {
          combinedProbs[, newClassName] <- rowSums(seuratObj@meta.data[, constituentProbCols, drop = FALSE])
        } else if (probabilityAggregation == "max") {
          combinedProbs[, newClassName] <- apply(seuratObj@meta.data[, constituentProbCols, drop = FALSE], 1, max)
        } else if (probabilityAggregation == "mean") {
          combinedProbs[, newClassName] <- rowMeans(seuratObj@meta.data[, constituentProbCols, drop = FALSE])
        }
      }
    }
    
    #add combined probabilities to metadata
    combinedProbsDf <- as.data.frame(combinedProbs)
    colnames(combinedProbsDf) <- paste0(modelName, "_Combined_prob_", colnames(combinedProbsDf), "_", modelVersion)
    seuratObj <- Seurat::AddMetaData(seuratObj, metadata = combinedProbsDf)
    
    #if relabelOrRecall == "recall", assign classes using combined probabilities
    if (relabelOrRecall == "recall") {
      probMat <- as.matrix(combinedProbsDf)
      if (recallMethod == "max") {
        maxIdx <- apply(probMat, 1, which.max)
        combinedClasses <- colnames(probMat)[maxIdx]
      } else if (recallMethod == "threshold") {
        # Assign the first class that exceeds threshold; if none, label Uncalled
        combinedClasses <- rep("Uncalled", nrow(probMat))
        for (i in seq_len(nrow(probMat))) {
          exceeds <- which(probMat[i, ] >= recallProbabilityThreshold)
          if (length(exceeds) > 0) {
            combinedClasses[i] <- colnames(probMat)[exceeds[1]]
          }
        }
      }
    }
  }
  
  #add combined classes to metadata
  #ensure the column name is tied to the model and version
  combinedClassDf <- data.frame(
    combinedClasses,
    row.names = colnames(seuratObj),
    stringsAsFactors = FALSE
  )
  colnames(combinedClassDf) <- outputFieldName
  seuratObj <- Seurat::AddMetaData(seuratObj, metadata = combinedClassDf)
  
  #print summary
  print(paste0("Combined ", length(unique(seuratObj@meta.data[[classFieldName]])), 
               " original classes into ", length(classMapping), " combined classes"))
  print(table(seuratObj@meta.data[[outputFieldName]]))
  
  #plot if reductions are available
  if (length(names(seuratObj@reductions)) > 0) {
    print(Seurat::DimPlot(seuratObj, group.by = outputFieldName))
  }
  
  return(seuratObj)
}


# ----------------------------------------------------------------------------
# Activation Class Mapping Registry
# ----------------------------------------------------------------------------

# initialize registry in package environment
if (is.null(pkg.env$ACTIVATION_CLASS_MAPPINGS)) {
  pkg.env$ACTIVATION_CLASS_MAPPINGS <- list()
}

.RegisterActivationClassMapping <- function(name, classMapping) {
  if (missing(name) || is.null(name) || nchar(name) == 0) {
    stop('Name must be a non-empty string')
  }
  if (name %in% names(pkg.env$ACTIVATION_CLASS_MAPPINGS)) {
    stop(paste0('Activation class mapping already registered: ', name))
  }

  # basic validation similar to CombineTcellActivationClasses
  if (!is.list(classMapping) || is.null(names(classMapping))) {
    stop('classMapping must be a named list')
  }
  if (any(names(classMapping) == '' | is.na(names(classMapping)))) {
    stop('All elements in classMapping must have non-empty names')
  }
  for (newClassName in names(classMapping)) {
    vals <- classMapping[[newClassName]]
    if (!is.character(vals)) {
      stop(paste0("Value for '", newClassName, "' must be a character vector"))
    }
    if (length(vals) == 0) {
      stop(paste0("Value for '", newClassName, "' cannot be an empty vector"))
    }
    if (any(is.na(vals) | vals == '')) {
      stop(paste0("Value for '", newClassName, "' contains empty or NA class names"))
    }
  }

  # duplicates across mappings are ambiguous
  allMapped <- unlist(classMapping, use.names = FALSE)
  dups <- allMapped[duplicated(allMapped)]
  if (length(dups) > 0) {
    stop(paste0('The following original classes are duplicated across mappings: ', paste0(unique(dups), collapse = ', ')))
  }

  pkg.env$ACTIVATION_CLASS_MAPPINGS[[name]] <- classMapping
}

#' @title GetActivationClassMapping
#'
#' @description Fetch a predefined T cell activation class mapping by key. These mappings
#' collapse model-specific fine-grained classes into broader categories and are suitable
#' inputs for CombineTcellActivationClasses().
#' @param name The key of the registered class mapping
#' @return A named list mapping new combined classes to character vectors of original classes.
#' Returns NULL and warns if the key is unknown.
#' @export
GetActivationClassMapping <- function(name) {
  if (!(name %in% names(pkg.env$ACTIVATION_CLASS_MAPPINGS))) {
    warning(paste0('Unknown activation class mapping: ', name))
    return(NULL)
  }
  return(pkg.env$ACTIVATION_CLASS_MAPPINGS[[name]])
}

# pre-register common mappings
.RegisterActivationClassMapping(
  'TcellActivation.Basic',
  list(
    'Th1-Tc1-like' = c('Th1_MIP1B.neg_CD137.neg', 'Th1_MIP1B.neg_CD137.pos', 'Th1_MIP1B.pos'),
    'Th2-Th17-like' = c('Th17'),
    'Bystander' = c('NonSpecificActivated_L1', 'NonSpecificActivated_L2'),
    'Bulk Tissue T cell' = c('Uncultured'),
    'Cultured T cell' = c('Cultured_Bystander_NoBFA', 'Cultured_Bystander_BFA')
  )
)

.RegisterActivationClassMapping(
  'TcellActivation.NomenclatureV1',
  list(
    'CD4posOrCD8pos_Th1_U_B_A_t_star' = c('Th1_MIP1B.neg_CD137.neg', 'Th1_MIP1B.neg_CD137.pos', 'Th1_MIP1B.pos'),
    'CD4pos_Th17_U_B_A_t_star' = c('Th17'),
    'CD4posOrCD8pos_Th1_U_B_A_t_O' = c('NonSpecificActivated_L1', 'NonSpecificActivated_L2'),
    'CD4posOrCD8pos_T_U_R_N_t_O' = c('Uncultured'),
    'CD4posOrCD8pos_T_U_B_N_t_O' = c('Cultured_Bystander_NoBFA', 'Cultured_Bystander_BFA')
  )
)
