# Predicts T cell activation using sPLS derived components and a trained logistic model on transformed variates

Predicts T cell activation using a trained model

## Usage

``` r
PredictTcellActivation(
  seuratObj,
  modelName = "GeneralizedTCR",
  combineClasses = TRUE
)
```

## Arguments

- seuratObj:

  The Seurat Object to be updated

- modelName:

  The model to use, either GeneralizedTCR, CD4, or CD8

- combineClasses:

  Most models were trained using finer grained classes than is generally
  needed. The default (combinedClasses=TRUE), will merge some similar
  classes. Set this to false to turn off this behavior.

## Value

A Seurat object with the sPLSDA scores and predicted probabilities added
to the metadata
