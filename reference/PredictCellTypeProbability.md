# Applies trained models to get celltype probabilities

Applies trained models to get celltype probabilities.

## Usage

``` r
PredictCellTypeProbability(
  seuratObj,
  models,
  fieldName = "RIRA_Consensus",
  batchSize = 20000,
  assayName = "RNA",
  minimum_probability = 0.5,
  minimum_delta = 0.25
)
```

## Arguments

- seuratObj:

  The Seurat Object to be updated

- models:

  A named vector of models, where the values are the modelName (for
  built-in model), or filePath to an RDS file. The names of the vector
  should be the cell-type label for cells scored positive by the
  classifier.

- fieldName:

  The name of the metadata column to store the result

- batchSize:

  To conserve memory, data will be chunked into batches of at most this
  many cells

- assayName:

  The assay holding gene expression data

- minimum_probability:

  The minimum probability for a confident cell type assignment

- minimum_delta:

  The minimum difference in probabilities necessary to call one celltype
  over another.
