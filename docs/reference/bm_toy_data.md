# Toy 4-country, 3-sector IO table for bmGVC

A small multi-country input–output data set used in bmGVC examples and
vignettes. It contains four countries (China, India, Japan, ROW) and
three sectors (Primary, Manufacturing, Service).

## Format

- bm_toy_Z:

  numeric matrix `12 x 12`

- bm_toy_Y:

  numeric matrix `12 x 4`

- bm_toy_VA:

  numeric vector of length 12

- bm_toy_X:

  numeric vector of length 12

- bm_toy_countries:

  character vector of length 4

- bm_toy_sectors:

  character vector of length 3

## Details

The data are stored in six objects:

- `bm_toy_Z`: 12 x 12 intermediate demand matrix

- `bm_toy_Y`: 12 x 4 final demand matrix

- `bm_toy_VA`: length-12 value-added vector

- `bm_toy_X`: length-12 gross output vector

- `bm_toy_countries`: character vector of length 4

- `bm_toy_sectors`: character vector of length 3

The ordering of industries is (China P,M,S; India P,M,S; Japan P,M,S;
ROW P,M,S).
