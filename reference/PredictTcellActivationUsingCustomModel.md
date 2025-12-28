# Predicts T cell activation using sPLS derived components and a trained logistic model on transformed variates

Predicts T cell activation using a trained model

## Usage

``` r
PredictTcellActivationUsingCustomModel(seuratObj, modelName, modelFile)
```

## Arguments

- seuratObj:

  The Seurat Object to be updated

- modelName:

  A name to use for this model. This will be the prefix of output
  columns

- modelFile:

  The path to an RDS with this model.

## Value

A Seurat object with the sPLSDA scores and predicted probabilities added
to the metadata
