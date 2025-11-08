# Applies a trained binary classifier to get per-cell probabilities.

Applies a trained model to get per-cell probabilities.

## Usage

``` r
ScoreCellsWithSavedModel(
  seuratObj,
  model,
  fieldToClass,
  batchSize = 20000,
  assayName = "RNA"
)
```

## Arguments

- seuratObj:

  The Seurat Object to be updated

- model:

  Either the full filepath to a model RDS file, or the name of a
  built-in model.

- fieldToClass:

  A list mapping the target field name in the seurat object to the
  classifier level. The latter is either numeric, or the string label.
  For example: list('CD4_T' = 1, 'CD8_T' = 2))

- batchSize:

  To conserve memory, data will be chunked into batches of at most this
  many cells

- assayName:

  The assay holding gene expression data
