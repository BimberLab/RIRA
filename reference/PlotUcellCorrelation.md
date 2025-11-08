# PlotUcellCorrelation

This will plot the correlation between a UCell score and a gene set of
interest

## Usage

``` r
PlotUcellCorrelation(seuratObj, toCalculate, assayName = "RNA")
```

## Arguments

- seuratObj:

  The seurat object

- toCalculate:

  A named list where each item is a character vector of genes

- assayName:

  The assay to use

## Value

A list with moduleName and the spearman correlation matrix
