# Assigns celltype label based on probabilities in the metadata

Assigns celltype label based on probabilities in the metadata

## Usage

``` r
AssignCellType(
  seuratObj,
  probabilityColumns,
  fieldName = "RIRA_Consensus",
  minimum_probability = 0.5,
  minimum_delta = 0.25
)
```

## Arguments

- seuratObj:

  The Seurat Object to be updated

- probabilityColumns:

  The set of columns containing probabilities for each classifier to
  include

- fieldName:

  The name of the metadata column to store the result

- minimum_probability:

  The minimum probability for a confident cell type assignment

- minimum_delta:

  The minimum difference in probabilities necessary to call one celltype
  over another.
