# Interprets the feature importance of each model

Interprets the feature importance of each model using DALEX and feature
importance permuation

## Usage

``` r
InterpretModels(output_dir = "./classifiers", plot_type = "ratio")
```

## Arguments

- output_dir:

  The output directory that TrainAllModels saved training data and
  models into

- plot_type:

  Argument to pass to model_parts(). Ratio or difference is recommended
  for large feature sets where the base model's AUC loss will be outside
  the plotting range.
