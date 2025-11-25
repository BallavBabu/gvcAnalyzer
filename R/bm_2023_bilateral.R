# R/bm_2023_bilateral.R

#' BM_2023 source-based bilateral decomposition of e_sr
#' @export
bm_2023_bilateral_source <- function(io, s, r) {
  stopifnot(inherits(io, "bm_io"))
  s <- bm_country_id(io, s)
  r <- bm_country_id(io, r)
  if (s == r) stop("s and r must be different")

  # TODO: implement DVAsource_sr, DDCsource_sr, FVAsource_sr, FDCsource_sr
  # using /s notation exactly as in your BM_2023 formulas.

  data.frame(
    exporter = io$countries[s],
    importer = io$countries[r],
    DVAsource_sr = NA_real_,
    DDCsource_sr = NA_real_,
    FVAsource_sr = NA_real_,
    FDCsource_sr = NA_real_
  )
}
