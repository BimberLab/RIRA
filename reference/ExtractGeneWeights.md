# ExtractGeneWeights

Extracts the loadings from a given component in a seurat object and
returns a data.frame suitable for ScoreUsingSavedComponent

## Usage

``` r
ExtractGeneWeights(
  seuratObj,
  componentNum,
  topNGenes = 100,
  bottomNGenes = 100,
  reduction = "pca"
)
```

## Arguments

- seuratObj:

  The seurat object

- componentNum:

  The component number to use

- topNGenes:

  The number of top genes to use (ranked by loading)

- bottomNGenes:

  The number of bottom genes to use (ranked by loading)

- reduction:

  The reduction to use
