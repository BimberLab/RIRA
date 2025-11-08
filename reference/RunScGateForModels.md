# Run scGate for models

Helper function to run scGate, iterating the provided models and
generating a consensus field

## Usage

``` r
RunScGateForModels(
  seuratObj,
  modelNames,
  min.cells = 30,
  assay = "RNA",
  pos.thr = 0.13,
  neg.thr = 0.13,
  ncores = 1,
  genes.blacklist = "default",
  labelRename = NULL,
  dropAmbiguousConsensusValues = FALSE,
  consensusModels = NULL
)
```

## Arguments

- seuratObj:

  The seurat object

- modelNames:

  A vector of model names to run. They are assumed to be non-overlapping
  populations

- min.cells:

  Passed directly to scGate::scGate. Stop iterating if fewer than this
  number of cells is left

- assay:

  Passed directly to scGate::scGate. Seurat assay to use

- pos.thr:

  Passed directly to scGate::scGate. Minimum UCell score value for
  positive signatures

- neg.thr:

  Passed directly to scGate::scGate. Maximum UCell score value for
  negative signatures

- ncores:

  Passed directly to scGate::scGate. Number of processors for parallel
  processing (requires future.apply)

- genes.blacklist:

  Passed directly to scGate::scGate. Genes blacklisted from variable
  features. The default loads the list of genes in
  scGate::genes.blacklist.default; you may deactivate blacklisting by
  setting genes.blacklist=NULL

- labelRename:

  An optional list that maps the model name to the final label that
  should be used in the seurat object. for exmaple: list(Tcell = 'T_NK',
  NK = 'T_NK'), would re-label cells classified as either 'Tcell' or
  'NK' by those models to one common label of T_NK

- dropAmbiguousConsensusValues:

  If true, any consensus calls that are ambiguous will be set to NA

- consensusModels:

  An optional list of model names to consider for the consensus call.
  This allows many models to be run, yet only consider a subset when
  creating the consensus call. This might be useful if some models
  overlap or produce false-positives.
