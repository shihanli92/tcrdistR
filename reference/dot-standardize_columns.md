# Rename columns from source format to tcrdistR canonical names

Rename columns from source format to tcrdistR canonical names

## Usage

``` r
.standardize_columns(df, col_map)
```

## Arguments

- df:

  A `data.frame`.

- col_map:

  A named list where names are target (canonical) column names and
  values are source column names present in `df`.

## Value

The `data.frame` with matching columns renamed.
