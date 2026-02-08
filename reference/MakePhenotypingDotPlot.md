# MakePhenotypingDotPlot

Creates a DotPlot using custom gene sets and attempts to coarsely group
gene sets by cell type.

## Usage

``` r
MakePhenotypingDotPlot(
  seuratObj,
  yField = "ClusterNames_0.2",
  scaled = T,
  gene_lists = c("Cytotoxicity.1", "EffectorCytokines"),
  assay = "RNA",
  scale.by = "size"
)
```

## Arguments

- seuratObj:

  A Seurat Object storing the count matrix to be used for phenotyping.

- yField:

  The grouping variable used to calculate the average expression of
  genes and the y axis of the DotPlot.

- scaled:

  A boolean defining whether to color dots by scaled expression or
  unscaled expression.

- gene_lists:

  A vector of gene lists (defined by .RegisterGeneSet) to be queried and
  their genes be plotted.

- assay:

  Which assay to use in the input seuratObj

- scale.by:

  Allow different scaling methods for dot size. 'radius' will
  de-emphasize lower/intermediately percent expressed genes.
