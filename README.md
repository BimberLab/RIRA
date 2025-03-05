[![R Build and Checks](https://github.com/bimberlab/RIRA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bimberlab/RIRA/actions/workflows/R-CMD-check.yaml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14976335.svg)](https://doi.org/10.5281/zenodo.14976335)

## Rhesus Immune Reference Atlas (RIRA): A multi-tissue single-cell landscape of immune cells

## Table of Contents
* [Overview](#overview)
* [Installation](#installation)
* [Download RIRA Data](https://bimberlab.github.io/RIRA/articles/GetRira.html)


### <a name = "overview">Overview</a>

RIRA is a multi-tissue single-cell RNA-seq atlas from Rhesus macaque cells, including data from 7 tissues and 47 animals. 
It provides highly annotated reference data, and a number of pipelines to assign cell type to your data using RIRA's built-in models.
Finally, it includes pipelines to replicate the process used to generate RIRA's reference labels, which might be useful for generating a similar reference for another species.      

The raw expression data are available, including a seurat object, through the NIH GEO database under accession [GSE277821](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277821).

[Please see our documentation for more detail and usage examples](https://bimberlab.github.io/RIRA/)


### <a name="installation">Installation</a>
```
# Make sure to update your Rprofile to include Bioconductor repos, such as adding this line to ~/.Rprofile:
local({options(repos = BiocManager::repositories())})

#Latest version:
devtools::install_github(repo = 'bimberlab/RIRA', ref = 'master', dependencies = TRUE, upgrade = 'always')
```

Pre-packaged Docker images with all needed dependencies installed can be found on our [GitHub Packages page](https://github.com/orgs/bimberlab/RIRA/pkgs/container/rira):

```
docker pull ghcr.io/bimberlab/rira:latest
```
