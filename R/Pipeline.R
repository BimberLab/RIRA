library(mlr3verse)
library(DALEX)
library(DALEXtra)
get_info_from_scGate <- function(model){
  info <- unique(no_antpres[, c("name", "signature")])
  celltypes <- info[, "name"]
  genes <- c()
  for (rowgenelist in info$signature){
    splitrowgenelist <- strsplit(rowgenelist, ";")[[1]]
    splitrowgenelist <- gsub("-$","", splitrowgenelist)
    genes <- c(genes, splitrowgenelist)
  }
  genes <- unique(genes)
  return(list(celltypes=celltypes, genes=genes))
}

create_training_matrix <- function(seuratObj, label_col, genelist){
  training_matrix <- as.data.frame(t(as.matrix(seuratObj@assays$RNA[genelist,])))
  training_matrix$label <- seuratObj@meta.data[,label_col]
  labels <- unique(training_matrix$label)
  numlabels <- length(labels) + 1
  for (celltype in labels){
    training_matrix[,celltype] <- ifelse(training_matrix$label==celltype,1,0)
  }
  numlabels
  return(list(training_matrix=training_matrix, numlabels=numlabels))
}

train_model <- function(training_matrix, label, numlabels, seedval=0){
  set.seed(seedval)
  markers_exp <- training_matrix
  colnames(markers_exp) <- gsub("-", ".", names(markers_exp))
  classification.data <- markers_exp[,1:(ncol(markers_exp)-numlabels)]
  genelist <- colnames(classification.data)
  labels_mat <- markers_exp[,(ncol(markers_exp)-numlabels+1):ncol(markers_exp)]
  
  classification.data[,label] <- as.factor(labels_mat[,label])
  colnames(classification.data) <- make.names(colnames(classification.data),unique = T)
  
  task <- TaskClassif$new(classification.data, id = "classifier", target = label, positive="1")
  learner = lrn("classif.ranger", importance = "permutation", num.trees=500, mtry=length(genelist), predict_type = "prob")
  set_threads(learner)
  train_set = sample(task$nrow, 0.8 * task$nrow)
  test_set = setdiff(seq_len(task$nrow), train_set)
  dalex_data <- markers_exp[train_set,1:(ncol(markers_exp)-numlabels)]
  dalex_y <- as.numeric(classification.data[train_set, label])
  
  pred = learner$train(task, row_ids = train_set)$predict(task, row_ids=test_set)
  confusion <- caret::confusionMatrix(factor(pred$response), factor(pred$truth))
  print(confusion$table)
  print(confusion$overall["Accuracy"])
  return(list(learner= learner, pred=pred, dalex_data=dalex_data, dalex_y=dalex_y))
}

train_all_models <- function(training_matrix, numlabels, seedval=0, print_dalex=T){
  markers_exp <- training_matrix
  labels_mat <- markers_exp[,(ncol(markers_exp)-numlabels+2):ncol(markers_exp)]
  models <- list()
  pred <- list()
  if (print_dalex==T) {
    dalex_res <- list()
    for (label in colnames(labels_mat)){
      print(label)
      temp <- train_model(training_matrix, label, numlabels, seedval)
      models[[label]] <- temp$learner
      pred[[label]] <- temp$pred
      dalex_res[[label]] <- produce_Dalex_results(temp$learner, temp$dalex_data, temp$dalex_y, label)
    }
    return(list(models=models, pred=pred, dalex_res=dalex_res))
  } else {
    for (label in colnames(labels_mat)){
      print(label)
      temp <- train_model(training_matrix, label, numlabels, seedval)
      models[[label]] <- temp$learner
      pred[[label]] <- temp$pred
    }
  }
  return(list(models=models, pred=pred))
}

produce_Dalex_results <- function(learner, dalex_data, dalex_y, label) {
  ranger_exp = explain_mlr3(learner,
                            data     = dalex_data,
                            y        = dalex_y,
                            B        = 1000,
                            label    = label,
                            colorize = FALSE)
  
  parts = model_parts(ranger_exp)
  print(plot(parts, max_vars=20, show_boxplots = FALSE))
  
  return(list(parts= parts, ranger_exp=ranger_exp))
}

genemat_from_seurat <- function(seuratObj, genelist){
  pred_matrix <- as.data.frame(t(as.matrix(seuratObj@assays$RNA[genelist,])))
  colnames(pred_matrix) <- gsub("-", ".", names(pred_matrix))
  return(pred_matrix)
}

make_predictions <- function(models, pred_matrix){
  celltypes <- names(model_res$models)
  df <- NULL
  for (celltype in celltypes){
    print(celltype)
    predvec <- model_res$models[[celltype]]$predict_newdata(pred_matrix)$prob[,1]
    if (is.null(df)) {
      df <- data.frame(predvec)
      names(df) <- celltype
    }
    else {
      df[[celltype]] <- predvec
    }
  }
  rownames(df) <- rownames(pred_matrix)
  return(df)
}

get_Top2 <- function(probdata) {
  df <- data.frame(probdata, check.names = FALSE)
  top <- c()
  toplabel <- c()
  secondlabel <- c()
  topval <- c()
  second <- c()
  for (i in 1:nrow(df)) {
    top[i] <- which.max(df[i, ])
    toplabel[i] <- colnames(df)[top[i]]
    topval[i] <- df[i, top[i]]
    df[i,top[i]] <- -1000
    second[i] <- max(df[i, ])
    secondlabel[i] <- colnames(df)[which.max(df[i, ])]
  }
  outdf <- data.frame("CellID" = rownames(df), check.names = FALSE)
  outdf$TopLabel <- toplabel
  outdf$Label2 <- secondlabel
  outdf$Highest <- topval
  outdf$Second <- second
  return(outdf)
}

assign_labels <- function(top2, minprob, mindelta){
  df <- data.frame(top2)
  df$Final <- ifelse(df$Highest > minprob & df$Highest-df$Second > mindelta, df$TopLabel, "Unknown")
  return(df)
}
