# Run Celltypist

Runs celltypist on a seurat object and stores the calls as metadata

## Usage

``` r
RunCellTypist(
  seuratObj,
  modelName = "Immune_All_Low.pkl",
  pThreshold = 0.5,
  minProp = 0,
  useMajorityVoting = TRUE,
  mode = "prob_match",
  extraArgs = c("--mode", mode, "--p-thres", pThreshold, "--min-prop", minProp),
  assayName = Seurat::DefaultAssay(seuratObj),
  columnPrefix = NULL,
  maxAllowableClasses = 6,
  minFractionToInclude = 0.01,
  minCellsToRun = 200,
  maxBatchSize = 1e+05,
  retainProbabilityMatrix = FALSE,
  runCelltypistUpdate = TRUE
)
```

## Arguments

- seuratObj:

  The seurat object

- modelName:

  The model name or path to celltypist model

- pThreshold:

  By default, this would be passed to the –p-thres argument. However, if
  you also provide extraArgs, this is ignored.

- minProp:

  By default, this would be passed to the –min-prop argument. However,
  if you also provide extraArgs, this is ignored.

- useMajorityVoting:

  If true, the celltypist –majority-voting option will be added.

- mode:

  Passed directly to the celltypist –mode argument. Either best_match or
  prob_match

- extraArgs:

  An optional list of additional arguments passed directly on the
  command line to cell typist

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

- runCelltypistUpdate:

  If true, –update-models will be run for celltypist prior to scoring
  cells.
