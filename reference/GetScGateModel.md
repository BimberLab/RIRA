# GetScGateModel

Returns the selected scGate model

## Usage

``` r
GetScGateModel(modelName, allowSCGateDB = TRUE)
```

## Arguments

- modelName:

  The name of the gate to return. See GetAvailableScGates() for a list
  of known gates

- allowSCGateDB:

  If true, this will search local models and the models provided by
  scGate::get_scGateDB()
