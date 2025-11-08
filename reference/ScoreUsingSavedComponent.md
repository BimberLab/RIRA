# ScoreUsingSavedComponent

Scores a seurat object using the gene weights from a saved component

## Usage

``` r
ScoreUsingSavedComponent(
  seuratObj,
  componentOrName,
  fieldName,
  assayName = "RNA",
  layer = "data"
)
```

## Arguments

- seuratObj:

  The seurat object

- componentOrName:

  Either a data frame with the columns 'feature' and 'weight', or the
  name of a previously saved component available in RIRA

- fieldName:

  The name of the field in which to save the resulting scores

- assayName:

  The assay to use

- layer:

  The layer to use
