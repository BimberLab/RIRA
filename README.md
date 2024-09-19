[![R Build and Checks](https://github.com/bimberlab/RIRA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bimberlab/RIRA/actions/workflows/R-CMD-check.yaml)

## Rhesus Immunome Reference Atlas (RIRA): A multi-tissue single-cell landscape of immune cells

## Table of Contents
* [Overview](#overview)
* [Example Usage](#usage)
* [Installation](#installation)


### <a name = "overview">Overview</a>

RIRA is a multi-tissue single-cell RNA-seq atlas from Rhesus macaque cells, including data from 7 tissues and 47 animals. 
It provides highly annotated reference data, and a number of pipelines to assign cell type to your data using RIRA's built-in models.
Finally, it includes pipelines to replicate the process used to generate RIRA's reference labels, which might be useful for generating a similar reference for another species.      

### <a name = "usage">Data</a>

RIRA's data was generated and labeled using the following process:
1) Merge data from all animals.
2) Perform strict QC and filtration, based on metrics such as per-cell RNA saturation
3) For initial phenotyping, we used scGate with a set of rhesus-adapted gates to divide cells. This initial pass is lossy, and will not label all cells.
4) Cells labeled by scGate were downsampled to equalize cells per class, and then fed as training data to celltypist. The resulting model (RIRA_Immune_v2) is available through this R package.
5) The above model will classify cells into course bins (T/NK, Bcell, Myeloid, and Other). For many analyses we find it useful to subdivide a dataset into these categories and then perform dimensionality reduction on these more homogenous subsets.
6) RIRA provides a handful of models trained to identify more specific subsets (such as differentiating CD4 vs CD8 T cells)

### <a name = "usage">Usage Examples</a>

```
# Use the built-in celltypist model to score cells according to course phenotypes (T/NK, Bcell, Myeloid, Other):
seuratObj <- RIRA::RunCellTypist(seuratObj, modelName = ‘RIRA_Immune_v2’)
Seurat::DimPlot(seuratObj, group.by = 'RIRA_Immune_v2.cellclass')

# Also the built-in celltypist model for finer scope T/NK or myeloid cells:
seuratObj <- RIRA::Classify_TNK(seuratObj)
Seurat::DimPlot(seuratObj, group.by = 'RIRA_TNK_v2.cellclass')

seuratObj <- RIRA::Classify_Myeloid(seuratObj)
Seurat::DimPlot(seuratObj, group.by = 'RIRA_Myeloid_v3.cellclass')

# USing sorted T cell data, we empirically defined a gene component that can classify T cells by naive-to-memory differentiation:
seuratObj <- ScoreUsingSavedComponent(seuratObj, componentOrName = 'Tcell_EffectorDifferentiation', fieldName = 'EDS')
Seurat::FeaturePlot(seuratObj, features = 'EDS')

```


### <a name="installation">Installation</a>

```{r}
# Make sure to update your Rprofile to include Bioconductor repos, such as adding this line to ~/.Rprofile:
local({options(repos = BiocManager::repositories())})

#Latest version:
devtools::install_github(repo = 'bimberlab/RIRA', ref = 'master', dependencies = TRUE, upgrade = 'always')
```

Pre-packaged Docker images with all needed dependencies installed can be found on our [GitHub Packages page](https://github.com/orgs/bimberlab/RIRA/pkgs/container/rira). We recommend using a specific release, which you can do using tags:

```
docker pull ghcr.io/bimberlab/rira:latest
```
