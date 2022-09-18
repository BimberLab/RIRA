[![R Build and Checks](https://github.com/bimberlabinternal/RIRA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bimberlabinternal/RIRA/actions/workflows/R-CMD-check.yaml)

## Rhesus Immunome Reference Atlas (RIRA): A multi-tissue single-cell landscape of immune cells

## Table of Contents
* [Overview](#overview)
* [Example Usage](#usage)
* [Installation](#installation)


### <a name = "overview">Overview</a>

RIRA is a multi-tissue single-cell RNAseq atlas from Rhesus macaque cells. 
It provides highly annotated reference data, and a number of pipelines to assign cell type to your data using RIRA's built-in models.
Finally, it includes pipelines to replicate the process used to generate RIRA's reference labels, which might be useful for generating a similar reference for another species.      

### <a name = "usage">Data</a>

RIRA's data was generated and labeled using the following process:
1) Merge data from XX 10x Genomics datasets, representing 8 tissues and XX rhesus macaques
2) Perform strict QC and filtration, based on metrics such as per-cell RNA saturation
3) For initial phenotyping, we used scGate with a set of rhesus-adapted gates to divide cells. This initial pass is lossy, and will not label all cells.
4) Cells labeled by scGate were downsampled to equalize cells per class, and then fed as training data to celltypist. The resulting model (RIRA_Immune_v1) is available through this R package.
5) The above model will classify cells into course bins (T/NK, Bcell, MoMacDC, and Other). For many analyses we find it useful to subdivide a dataset into these categories and then perform dimensionality reduction on these more homogenous subsets.
6) RIRA provides a handful of models trained to identify more specific subsets (such as differentiating CD4 vs CD8 T cells)

### <a name = "usage">Usage Examples</a>

```
# Use the built-in celltypist model to score cells according to course phenotypes (T/NK, Bcell, MoMacDC, Other):
seuratObj <- RIRA::RunCellTypist(seuratObj, modelName = ‘RIRA_Immune_v1’)

# Use the built-in celltypist model for course phenotypes, followed by higher-resolution labeling of T/NK cells:
seuratObj <- RIRA::ClassifyCells(seuratObj, primaryModel = ‘RIRA_Immune_v1’, subsetModels = list(
    T_NK = 'RIRA_CD4vCD8_LR'
)

```


### <a name="installation">Installation</a>

```{r}
# Make sure to update your Rprofile to include Bioconductor repos, such as adding this line to ~/.Rprofile:
local({options(repos = BiocManager::repositories())})

#Latest version:
devtools::install_github(repo = 'bimberlabinternal/RIRA', ref = 'master', dependencies = TRUE, upgrade = 'always')
```

Pre-packaged Docker images with all needed dependencies installed can be found on our [GitHub Packages page](https://github.com/orgs/BimberLabInternal/RIRA/pkgs/container/rira). We recommend using a specific release, which you can do using tags:

```
docker pull ghcr.io/bimberlabinternal/rira:latest
```
