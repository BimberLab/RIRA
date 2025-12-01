# GetActivationClassMapping

Fetch a predefined T cell activation class mapping by key. These
mappings collapse model-specific fine-grained classes into broader
categories and are suitable inputs for CombineTcellActivationClasses().

## Usage

``` r
GetActivationClassMapping(name)
```

## Arguments

- name:

  The key of the registered class mapping

## Value

A named list mapping new combined classes to character vectors of
original classes. Returns NULL and warns if the key is unknown.
