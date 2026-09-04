# RecoverUnassignedCells

This step recovers unassigned cells by inspecting each cluster and
assigning Unassigned cells to the match the cluster majority

## Usage

``` r
RecoverUnassignedCells(
  seuratObj,
  classField = "RIRA_Immune_v2.cellclass",
  groupField = "ClusterNames_0.2",
  targetField = "RIRA_Immune_v2.cellclass.recovered",
  unassignedValues = c("Unassigned", "Unknown"),
  minClusterProp = 0.6
)
```

## Arguments

- seuratObj:

  The seurat object

- classField:

  The name of the field holding the cell type call

- groupField:

  The field on which to group, such as the cluster ID

- targetField:

  The field that will store the result

- unassignedValues:

  A list of string values that denote Unassigned/Unknown cells

- minClusterProp:

  If at least this proportion of the group is one class, unassigned
  cells will be assigned as this class

## Value

The updated seurat object
