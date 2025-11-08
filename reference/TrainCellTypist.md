# Train Celltypist

Train celltypist using the data in a seurat object and saves the
resulting model

## Usage

``` r
TrainCellTypist(
  seuratObj,
  labelField,
  modelFile,
  minCellsPerClass = 20,
  assayName = Seurat::DefaultAssay(seuratObj),
  tempFileLocation = NULL,
  dropAmbiguousLabelValues = TRUE,
  excludedClasses = NULL,
  featureInclusionList = NULL,
  featureExclusionList = NULL
)
```

## Arguments

- seuratObj:

  The seurat object

- labelField:

  The field in seuratObj@meta.data holding the labels for training

- modelFile:

  The path to save the model

- minCellsPerClass:

  If provided, any classes (and corresponding cells) with fewer than
  this many cells will be dropped from the training data

- assayName:

  The name of the assay to use

- tempFileLocation:

  The location where temporary files (like the annData version of the
  seurat object), will be written.

- dropAmbiguousLabelValues:

  If true, and label value with a comma will be dropped.

- excludedClasses:

  A vector of labels to discard.

- featureInclusionList:

  If provided, the input count matrix will be subset to just these
  features. If used, Seurat::NormalizeData will be re-run.

- featureExclusionList:

  If provided, the input count matrix will be subset to remove these
  features. If used, Seurat::NormalizeData will be re-run.
