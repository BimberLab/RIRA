# Download RIRA data from GEO

This will download the RIRA raw data from GEO, as a Seurat object

## Usage

``` r
DownloadFromGEO(outfile = NULL, subset = c("All", "T_NK", "Myeloid", "Bcell"))
```

## Arguments

- outfile:

  This location to initially download the raw data. If omitted or NULL,
  it will download the RDS to the temporary file folder and delete it on
  completion. However, because the download is large it may be useful to
  download to a location you plan to keep.

- subset:

  The cell subset to download, either All, T_NK, Myeloid, BCell or Other
