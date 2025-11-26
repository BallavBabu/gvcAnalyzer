# Build a bm_io object from IO table blocks

Build a bm_io object from IO table blocks

## Usage

``` r
bm_build_io(Z, Y, VA, X, countries, sectors)
```

## Arguments

- Z:

  Intermediate demand matrix (GN x GN).

- Y:

  Final demand matrix. Can be (GN x G) OR (GN x (G \* FD_categories)).

- VA:

  Value added. Can be a vector (length GN) or matrix (Rows x GN).

- X:

  Output vector (length GN).

- countries:

  Character vector of country names/codes (length G).

- sectors:

  Character vector of sector names/codes (length N).

## Value

An object of class `"bm_io"`.
