[![R Build and Checks](https://github.com/bimberlabinternal/RIRA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bimberlabinternal/RIRA/actions/workflows/R-CMD-check.yaml)

## Rhesus Immunome Reference Atlas (RIRA): A multi-tissue single-cell landscape of immune cells

Code for utility of RIRA such as classification of cell types/phenotypes

## Table of Contents
* [Overview](#overview)
* [Example Usage](#usage)
* [Installation](#installation)


### <a name = "overview">Overview</a>

### <a name = "usage">Usage</a>

RIRA currently requires you to specify a local folder holding RIRA data. This is accomplished using:

```
SetAtlasDir('/path/to/RIRA')
```

Within this folder, there should be top-level folders corresponding to each version. Within each of these, it expects the following:

```
RIRA_Data/
    1.0.0/
        counts/
            barcodes.tsv.gz
            features.tsv.gz
            matrix.mtx.gz
        meta.csv
        cFIT/
            ref.rds
        
```

The counts folder is passed to DropletUtils::read10xCounts to create the raw count data. meta.csv is the meta.data data frame. The file ref.rds in ./cFIT is the integrated cFIT reference.

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
