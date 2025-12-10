# CalculateUCellScores

This will calculate a handful of standardized UCell scores for a seurat
object

## Usage

``` r
CalculateUCellScores(
  seuratObj,
  forceRecalculate = FALSE,
  seed = GetSeed(),
  ncores = 1,
  assayName = "RNA",
  storeRanks = TRUE,
  plotCor = TRUE,
  dropAllExistingUcells = FALSE
)
```

## Arguments

- seuratObj:

  The seurat object

- forceRecalculate:

  If true, UCell will always be re-run, even if the field is already
  present.

- seed:

  If non-null, set.seed() will be called prior to running UCell

- ncores:

  The number of core to use with UCell::AddModuleScore_UCell

- assayName:

  The assay from which to calculate UCell scores.

- storeRanks:

  Passed directly to UCell::AddModuleScore_UCell. Increases object size
  but makes future calculations quicker.

- plotCor:

  If true, a plot of correlations between the UCell score and each
  component gene will be shown

- dropAllExistingUcells:

  If true, any existing column ending in UCell will be removed
