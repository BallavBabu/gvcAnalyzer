# Introduction to gvcAnalyzer: BM 2023 and BM 2025 GVC Decomposition

## Introduction

The **gvcAnalyzer** package implements two complementary methodologies
for measuring global value chain (GVC) participation:

1.  **BM 2023**: Measuring what matters in value-added trade (Bilateral
    Focus).
2.  **BM 2025**: Economic consequences of trade and global value chain
    integration (Tripartite & Output Focus).

This vignette demonstrates both frameworks using a four-country,
three-sector input–output table. The example includes China, India,
Japan, and Rest of World (ROW), each with Primary, Manufacturing, and
Service sectors.

### References

Borin, A., & Mancini, M. (2023). Measuring what matters in value-added
trade. Economic Systems Research, 35(4), 586–613

Borin, A., Mancini, M., & Taglioni, D. (2025). Economic consequences of
trade and global value chain integration: A measurement perspective. The
World Bank Economic Review

## Data and Setup

### Toy Multi-Regional Input–Output Table

The package includes a minimal MRIO table for demonstration purposes:

``` r
# Countries and sectors
bm_toy_countries
#> [1] "China" "India" "Japan" "ROW"
bm_toy_sectors
#> [1] "Primary"       "Manufacturing" "Service"

# Dimensions
cat("Intermediate flows (Z):", paste(dim(bm_toy_Z), collapse = " × "), "\n")
#> Intermediate flows (Z): 12 × 12
cat("Final demand (Y):", paste(dim(bm_toy_Y), collapse = " × "), "\n")
#> Final demand (Y): 12 × 4
cat("Value added (VA):", length(bm_toy_VA), "industries\n")
#> Value added (VA): 12 industries
cat("Gross output (X):", length(bm_toy_X), "industries\n")
#> Gross output (X): 12 industries
```

### Building the IO Object

All **gvcAnalyzer** functions operate on a `bm_io` object that contains
the structured IO data and derived matrices:

``` r
io <- bm_build_io(
  Z         = bm_toy_Z,
  Y         = bm_toy_Y,
  VA        = bm_toy_VA,
  X         = bm_toy_X,
  countries = bm_toy_countries,
  sectors   = bm_toy_sectors
)

# Structure
io$G    # number of countries
#> [1] 4
io$N    # number of sectors
#> [1] 3
io$GN   # total industries
#> [1] 12
```

The `bm_io` object includes:

- Technical coefficient matrix **A**
- Global Leontief inverse **B**
- Value-added coefficient vector **v**
- Country-level domestic Leontief inverses **L**

## BM 2023: Detailed Bilateral Decomposition

The BM 2023 framework focuses on precise accounting of value-added in
bilateral trade, distinguishing between source and sink perspectives.

### Exporter-Level Decomposition

For each exporting country *s*, gross exports are decomposed into:

- **DVA_s**: Domestic value added
- **FVA_s**: Foreign value added
- **DDC_s**: Double-counted domestic content
- **FDC_s**: Double-counted foreign content

``` r
# Single country
bm_2023_exporter_total(io, 1) # Using index 1 for China
#>   country   DVA_s   DDC_s    FVA_s   FDC_s    EX_s
#> 1   China 2144685 8288.95 341856.2 1346.26 2496178

# All countries
bm_2023_exporter_total_all(io)
#>   country     DVA_s      DDC_s     FVA_s      FDC_s    EX_s
#> 1   China 2144685.3  8288.9498 341856.18 1346.25971 2496178
#> 2   India  424505.3   249.2700  72383.77   46.17244  497185
#> 3   Japan  660082.5   637.2357  94197.73   92.00585  755010
#> 4     ROW 2767872.8 13273.2642  94653.40  456.62833 2876256
```

### Bilateral Decomposition

#### Source-Based Perspective

The source-based decomposition tracks value added by its country of
origin in bilateral trade flows:

``` r
# Exports from China to India
bm_2023_bilateral_source(io, 1, 2)
#>   exporter importer DVAsource_sr DDCsource_sr FVAsource_sr FDCsource_sr EX_sr
#> 1    China    India     77542.28      300.136     12376.79     48.74709 90268
```

#### Pure Bilateral Flows

The pure decomposition isolates direct bilateral value-added flows,
excluding third-country effects:

``` r
bm_2023_bilateral_pure(io, 1, 2)
#>   exporter importer DVA_star_sr DDC_star_sr FVA_star_sr FDC_star_sr EX_sr
#> 1    China    India    77838.34    4.075798    12424.87   0.6656287 90268
```

## BM 2025: Tripartite & Output Decomposition

The BM 2025 framework introduces the “Tripartite” concept (Forward,
Backward, Two-Sided) and extends GVC measurement to the production side.

### Tripartite GVC Trade Decomposition

The tripartite decomposition classifies bilateral exports into
traditional trade and three types of GVC-related trade:

- **DAVAX**: Traditional trade (one-crossing domestic value added
  absorbed abroad)
- **GVC_PF**: Pure-forward GVC trade (domestic value added re-exported
  by the partner)
- **GVC_TS**: Two-sided GVC trade (simultaneous use of foreign inputs
  and re-export)
- **GVC_PB**: Pure-backward GVC trade (foreign value added in exports)

``` r
# Single bilateral pair (China -> India)
bm_2025_tripartite_trade(io, 1, 2)
#>   exporter importer  E_sr DAVAX_sr   GVC_sr   GVC_PF   GVC_TS   GVC_PB
#> 1    China    India 90268 69106.63 21161.37 8435.702 1403.456 11322.22
```

#### Aggregate Trade Components

Summing over all destination countries provides exporter-level GVC trade
components:

``` r
trade_comp <- bm_2025_trade_exporter(io)
trade_comp
#>   exporter     E_s     GVC_s  GVC_PF_s  GVC_TS_s  GVC_PB_s
#> 1    China 2496178 438493.93  87002.55 14403.632 337087.75
#> 2    India  497185  88754.79  16075.58  2975.461  69703.75
#> 3    Japan  755010 138716.78  43789.81  6550.941  88376.04
#> 4      ROW 2876256 547792.94 439409.60 17276.110  91107.23
```

#### GVC Participation Indicators

From these components, we compute trade-based participation measures:

``` r
trade_meas <- bm_2025_trade_measures(io)
trade_meas
#>   exporter     E_s     GVC_s  GVC_PF_s  GVC_TS_s  GVC_PB_s share_GVC_trade
#> 1    China 2496178 438493.93  87002.55 14403.632 337087.75       0.1756661
#> 2    India  497185  88754.79  16075.58  2975.461  69703.75       0.1785146
#> 3    Japan  755010 138716.78  43789.81  6550.941  88376.04       0.1837284
#> 4      ROW 2876256 547792.94 439409.60 17276.110  91107.23       0.1904535
#>   share_PF_trade share_TS_trade share_PB_trade forward_trade
#> 1      0.1984122     0.03284796      0.7687398    -0.5703276
#> 2      0.1811235     0.03352451      0.7853519    -0.6042284
#> 3      0.3156778     0.04722529      0.6370969    -0.3214191
#> 4      0.8021454     0.03153766      0.1663169     0.6358285
```

Key indicators:

- **share_GVC_trade**: Share of exports that are GVC-related
- **share_PF_trade**, **share_TS_trade**, **share_PB_trade**:
  Composition of GVC trade
- **forward_trade**: Forward orientation index = (GVC_PF − GVC_PB) / GVC

A positive `forward_trade` indicates upstream positioning (supplier of
intermediates), while negative values indicate downstream positioning
(assembler using foreign inputs).

## BM 2025: Output-Based GVC Decomposition

The BM 2025 framework measures GVC participation from the production
side, decomposing gross output rather than exports.

### Country-Level Output Components

For each country *s*, gross output is decomposed into:

- **DomX**: Purely domestic production (value added never crossing
  borders)
- **TradX**: Traditional one-crossing trade
- **GVC_PF_X**: Pure-forward GVC output
- **GVC_PB_X**: Pure-backward GVC output
- **GVC_TSImp**: Two-sided GVC output via imported intermediates
- **GVC_TSDom**: Two-sided GVC output via domestic intermediates
- **GVC_TS_X**: Total two-sided GVC output
- **GVC_X**: Total GVC-related output

``` r
out_comp <- bm_2025_output_components(io)
out_comp
#>   country  GVC_PF_X  GVC_PB_X GVC_TSImp GVC_TSDom  GVC_TS_X     GVC_X     DomX
#> 1   China  87001.08 226301.71 1984093.0 145326.53 2129419.5 2442722.3 12077920
#> 2   India  16075.05  46234.14  341713.9  13620.80  355334.7  417643.9  2151580
#> 3   Japan  43789.25  58926.79  527442.6  36361.99  563804.6  666520.6  4151500
#> 4     ROW 439409.35 254812.04 1823517.6 449083.39 2272601.0 2966822.4 57671452
#>      TradX   X_total
#> 1 22151856  36672498
#> 2  2270374   4839598
#> 3  3932253   8750273
#> 4 55587908 116226183
```

The identity holds: **X_total = DomX + TradX + GVC_X**

### Output-Based Participation Measures

From the output components, we derive participation indicators analogous
to the trade-based measures:

``` r
out_meas <- bm_2025_output_measures(io)
out_meas
#>   country  GVC_PF_X  GVC_PB_X GVC_TSImp GVC_TSDom  GVC_TS_X     GVC_X     DomX
#> 1   China  87001.08 226301.71 1984093.0 145326.53 2129419.5 2442722.3 12077920
#> 2   India  16075.05  46234.14  341713.9  13620.80  355334.7  417643.9  2151580
#> 3   Japan  43789.25  58926.79  527442.6  36361.99  563804.6  666520.6  4151500
#> 4     ROW 439409.35 254812.04 1823517.6 449083.39 2272601.0 2966822.4 57671452
#>      TradX   X_total share_GVC_output share_PF_output share_TS_output
#> 1 22151856  36672498       0.06660911      0.03561645       0.8717403
#> 2  2270374   4839598       0.08629724      0.03848983       0.8508079
#> 3  3932253   8750273       0.07617141      0.06569826       0.8458922
#> 4 55587908 116226183       0.02552628      0.14810774       0.7660051
#>   share_PB_output forward_output
#> 1      0.09264324    -0.05702680
#> 2      0.11070228    -0.07221245
#> 3      0.08840956    -0.02271129
#> 4      0.08588719     0.06222055
```

Key indicators:

- **share_GVC_output**: GVC-related output as a share of total output
- **share_PF_output**, **share_TS_output**, **share_PB_output**:
  Composition of GVC output
- **forward_output**: Output-based forward orientation index

### Sector-Level Decomposition

The BM 2025 framework extends to country–sector pairs by proportional
allocation based on sectoral output shares:

``` r
out_comp_sec <- bm_2025_output_components_sector(io)
head(out_comp_sec, 9)
#>   country        sector      X_i     DomX_i     TradX_i  GVC_PF_Xi   GVC_PB_Xi
#> 1   China       Primary  3281502 1104751.84  1886504.79 11065.2952   2197.2355
#> 2   China Manufacturing 16650390 3024032.37 11438638.98 44890.1216 158581.0612
#> 3   China       Service 16740606 7949135.51  8173835.06 31045.6670  32428.5510
#> 4   India       Primary   637990  432358.86   184911.01  3253.3344    528.8139
#> 5   India Manufacturing  1532351  293393.55  1011596.14  4467.9004  30579.2489
#> 6   India       Service  2669257 1425827.50  1138779.59  8353.8112  14137.8388
#> 7   Japan       Primary   122077   46272.56    65053.01   608.7586    267.9981
#> 8   Japan Manufacturing  2653472  721062.50  1599990.69 20695.8292  36667.2447
#> 9   Japan       Service  5974724 3384164.82  2368544.27 22484.6602  19152.4677
#>   GVC_TSImp_i GVC_TSDom_i   GVC_TS_Xi     GVC_Xi
#> 1  250445.495   26537.340  276982.835  290245.37
#> 2 1924357.616   59889.847 1984247.463 2187718.65
#> 3  495261.865   58899.345  554161.210  617635.43
#> 4   12386.388    4551.597   16937.984   20720.13
#> 5  189469.273    2844.887  192314.161  227361.31
#> 6   75933.949    6224.317   82158.267  104649.92
#> 7    8852.171    1022.506    9874.677   10751.43
#> 8  261587.673   13468.070  275055.743  332418.82
#> 9  158506.371   21871.418  180377.789  222014.92
```

## Summary

This vignette introduced the **gvcAnalyzer** package and demonstrated:

1.  Construction of the `bm_io` object from multi-regional input–output
    data
2.  BM 2023 bilateral decompositions for precise accounting
3.  BM 2025 tripartite trade decomposition for GVC roles
4.  BM 2025 output-based decomposition for production-side analysis

For empirical applications, users can replace the toy data with full
MRIO tables (e.g., WIOD, OECD-ICIO, ADB-MRIO) following the same
workflow.

## Session Information

``` r
sessionInfo()
#> R version 4.5.1 (2025-06-13)
#> Platform: aarch64-apple-darwin20
#> Running under: macOS Sequoia 15.6
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> time zone: Asia/Tokyo
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] gvcAnalyzer_0.1.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.38     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] Matrix_1.7-4      xfun_0.54         lattice_0.22-7    cachem_1.1.0     
#>  [9] knitr_1.50        htmltools_0.5.8.1 rmarkdown_2.30    lifecycle_1.0.4  
#> [13] cli_3.6.5         grid_4.5.1        sass_0.4.10       pkgdown_2.2.0    
#> [17] textshaping_1.0.4 jquerylib_0.1.4   systemfonts_1.3.1 compiler_4.5.1   
#> [21] rstudioapi_0.17.1 tools_4.5.1       ragg_1.5.0        bslib_0.9.0      
#> [25] evaluate_1.0.5    yaml_2.3.10       jsonlite_2.0.0    rlang_1.1.6      
#> [29] fs_1.6.6          htmlwidgets_1.6.4
```
