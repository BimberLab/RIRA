# PlotMarkerSeries

Iteratively plots a set of markers

## Usage

``` r
PlotMarkerSeries(
  seuratObj,
  features,
  reductions = c("umap"),
  title = NULL,
  setSize = 4
)
```

## Arguments

- seuratObj:

  The seurat object

- features:

  A vector of feature names

- reductions:

  The reductions to plot

- title:

  An optional title of this plot series

- setSize:

  The maximum number of features to include per FeaturePlot
