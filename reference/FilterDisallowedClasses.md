# Filter Disallowed Classes

This is used to flag cells with high UCell score combinations generally
indiciative of contamination or doublets

## Usage

``` r
FilterDisallowedClasses(
  seuratObj,
  sourceField = "RIRA_Immune_v2.majority_voting",
  outputFieldName = "DisallowedUCellCombinations",
  ucellCutoff = 0.2,
  disallowedClasses = list(T_NK = c("Bcell.RM_UCell", "Myeloid.RM_UCell",
    "Erythrocyte.RM_UCell", "Platelet.RM_UCell", "NeutrophilLineage.RM_UCell"), Myeloid =
    c("Bcell.RM_UCell", "Tcell.RM_UCell", "NK.RM_UCell", "Erythrocyte.RM_UCell",
    "Platelet.RM_UCell"), Bcell = c("Tcell.RM_UCell", "NK.RM_UCell", "Myeloid.RM_UCell",
    "Erythrocyte.RM_UCell", "Platelet.RM_UCell", "NeutrophilLineage.RM_UCell",
    "Complement.RM_UCell"))
)
```

## Arguments

- seuratObj:

  The seurat object

- sourceField:

  The name of the field on which to compare

- outputFieldName:

  The name of the field to store the results

- ucellCutoff:

  Any cells expressing the disallowed UCell above this value will be
  flagged

- disallowedClasses:

  This is a list where the names are the cell classes (which should
  match levels in sourceField), and values are a vector of UCell field
  names.
