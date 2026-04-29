# Validate and clean a TCR data.frame

Checks that required columns exist for the specified chain type, ensures
no NAs in critical columns, coerces factors to character, and returns
the cleaned data.frame.

## Usage

``` r
.validate_tcr_df(df, chains = "AB")
```

## Arguments

- df:

  A `data.frame` of TCR clonotypes.

- chains:

  Character string. One of `"AB"`, `"A"`, or `"B"`.

## Value

The cleaned `data.frame`.
