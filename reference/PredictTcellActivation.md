# Predicts T cell activation using sPLS derived components and a trained logistic model on transformed variates

Predicts T cell activation using a trained model

## Usage

``` r
PredictTcellActivation(seuratObj, model = NULL, modelList = NULL)
```

## Arguments

- seuratObj:

  The Seurat Object to be updated

- model:

  The trained sPLSDA model to use for prediction. This can be a file
  path to an RDS file, or a built-in model name.

- modelList:

  A list of trained sPLSDA models to use for prediction. This can be a
  list of file paths to RDS files, or built-in model names.

## Value

A Seurat object with the sPLSDA scores and predicted probabilities added
to the metadata
