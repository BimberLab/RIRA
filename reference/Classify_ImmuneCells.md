# Classify Bulk Immune cells

Runs celltypist using the RIRA bulk immune model to score cells using
CellTypist with optimized parameters.

## Usage

``` r
Classify_ImmuneCells(
  seuratObj,
  assayName = Seurat::DefaultAssay(seuratObj),
  columnPrefix = "RIRA_Immune_v2.",
  maxAllowableClasses = 6,
  minFractionToInclude = 0.01,
  minCellsToRun = 200,
  maxBatchSize = 6e+05,
  retainProbabilityMatrix = FALSE,
  filterDisallowedClasses = TRUE
)
```

## Arguments

- seuratObj:

  The seurat object

- assayName:

  The name of the assay to use. Others will be dropped

- columnPrefix:

  A prefix that will be added to the beginning of the resulting columns,
  added the seurat@meta.data

- maxAllowableClasses:

  Celltypist can assign a cell to many classes, creating extremely long
  labels. Any cell with more than this number of labels will be set to
  NA

- minFractionToInclude:

  If non-null, any labels with fewer than this fraction of cells will be
  set to NA.

- minCellsToRun:

  If the input seurat object has fewer than this many cells, NAs will be
  added for all expected columns and celltypist will not be run.

- maxBatchSize:

  If more than this many cells are in the object, it will be split into
  batches of this size and run in serial.

- retainProbabilityMatrix:

  If true, the celltypist probability_matrix with per-class
  probabilities will be stored in meta.data

- filterDisallowedClasses:

  If true, this will run FilterDisallowedClasses() on the output.
