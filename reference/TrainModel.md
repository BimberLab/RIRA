# Creates a binary classifier to classify cells within a Seurat object

Creates a binary classifier to classify cells

## Usage

``` r
TrainModel(
  training_matrix,
  celltype,
  hyperparameter_tuning = F,
  learner = "classif.ranger",
  inner_resampling = "cv",
  outer_resampling = "cv",
  inner_folds = 4,
  inner_ratio = 0.8,
  outer_folds = 3,
  outer_ratio = 0.8,
  n_models = 20,
  n_cores = NULL
)
```

## Arguments

- training_matrix:

  A matrix (counts or data layer) provided by TrainModelsFromSeurat

- celltype:

  The celltype (provided by TrainModelsFromSeurat) used as classifier's
  positive prediction

- hyperparameter_tuning:

  logical that determines whether or not hyperparameter tuning should be
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
