# RIRA

## Rhesus Immunome Reference Atlas (RIRA): A multi-tissue single-cell landscape of immune cells

Code for utility of RIRA such as classification of cell types/phenotypes


### Usage:

RIRA currently requires you to specify a local folder holding RIRA data. This is accomplished using:

```
SetAtlasDir('/path/to/RIRA)
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