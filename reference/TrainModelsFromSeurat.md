# Wrapper function for TrainModel to train a suite of binary classifiers for each cell type present in the data

Wrapper function for TrainModel to train a suite of binary classifiers
for each cell type present in the data

## Usage

``` r
TrainModelsFromSeurat(
  seuratObj,
  celltype_column,
  assay = "RNA",
  layer = "data",
  output_dir = "./classifiers",
  hyperparameter_tuning = F,
  learner = "classif.ranger",
  inner_resampling = "cv",
  outer_resampling = "cv",
  inner_folds = 4,
  inner_ratio = 0.8,
  outer_folds = 3,
  outer_ratio = 0.8,
  n_models = 20,
  n_cores = NULL,
  gene_list = NULL,
  gene_exclusion_list = NULL,
  verbose = TRUE,
  min_cells_per_class = 20
)
```

## Arguments

- seuratObj:

  The Seurat Object to be updated

- celltype_column:

  The metadata column containing the celltypes. One classifier will be
  created for each celltype present in this column.

- assay:

  SeuratObj assay containing the desired count matrix/metadata

- layer:

  Layer containing the count data. Should be restricted to counts, data,
  or scale.data.

- output_dir:

  The directory in which models, metrics, and training data will be
  saved.

- hyperparameter_tuning:

  Logical that determines whether or not hyperparameter tuning should be
  performed.

- learner:

  The mlr3 learner that should be used. Currently fixed to
  "classif.ranger" if hyperparameter tuning is FALSE. Otherwise,
  "classif.xgboost" and "classif.ranger" are supported.

- inner_resampling:

  The resampling strategy that is used for hyperparameter optimization.
  Holdout ("hout" or "holdout") and cross validation ("cv" or
  "cross-validation") are supported.

- outer_resampling:

  The resampling strategy that is used to determine overfitting. Holdout
  ("hout" or "holdout") and cross validation ("cv" or
  "cross-validation") are supported.

- inner_folds:

  The number of folds to be used for inner_resampling if
  cross-valdiation is performed.

- inner_ratio:

  The ratio of training to testing data to be used for inner_resampling
  if holdout resampling is performed.

- outer_folds:

  The number of folds to be used for outer_resampling if
  cross-valdiation is performed.

- outer_ratio:

  The ratio of training to testing data to be used for inner_resampling
  if holdout resampling is performed.

- n_models:

  The number of models to be trained during hyperparameter tuning. The
  model with the highest accuracy will be selected and returned.

- n_cores:

  If non-null, this number of workers will be used with future::plan

- gene_list:

  If non-null, the input count matrix will be subset to these features

- gene_exclusion_list:

  If non-null, the input count matrix will be subset to drop these
  features

- verbose:

  Whether or not to print the metrics data for each model after
  training.

- min_cells_per_class:

  If provided, any classes (and corresponding cells) with fewer than
  this many cells will be dropped from the training data
