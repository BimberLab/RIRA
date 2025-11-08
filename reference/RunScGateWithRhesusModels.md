# Run scGate using Rhesus macaque models

Helper function to run scGate, iterating the provided models and
generating a consensus field

## Usage

``` r
RunScGateWithRhesusModels(
  seuratObj,
  min.cells = 30,
  assay = "RNA",
  pos.thr = 0.13,
  neg.thr = 0.13,
  ncores = 1,
  genes.blacklist = "default",
  dropAmbiguousConsensusValues = FALSE
)
```

## Arguments

- seuratObj:

  The seurat object

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

- dropAmbiguousConsensusValues:

  If true, any consensus calls that are ambiguous will be set to NA
