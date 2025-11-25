# R/bm_2023_wrappers.R
#
# Convenience wrappers for BM_2023:
#   - All exporters (national totals)
#   - All bilateral pairs (source-based and pure /sr)

#' BM_2023 exporter totals for all countries
#'
#' Applies \code{bm_2023_exporter_total()} to every country in the IO object.
#'
#' @param io bm_io object
#'
#' @return data.frame with one row per exporter:
#'   country, DVA_s, DDC_s, FVA_s, FDC_s, EX_s
#' @export
bm_2023_exporter_total_all <- function(io) {
  stopifnot(inherits(io, "bm_io"))
  do.call(
    rbind,
    lapply(seq_len(io$G), function(s) bm_2023_exporter_total(io, s))
  )
}

#' BM_2023 source-based bilateral decomposition for all pairs
#'
#' Applies \code{bm_2023_bilateral_source()} to all ordered pairs (s,r), s ≠ r.
#'
#' @param io bm_io object
#'
#' @return data.frame with one row per exporter-importer pair:
#'   exporter, importer, DVAsource_sr, DDCsource_sr, FVAsource_sr, FDCsource_sr, EX_sr
#' @export
bm_2023_bilateral_source_all <- function(io) {
  stopifnot(inherits(io, "bm_io"))
  G <- io$G

  records <- list()
  k <- 1L

  for (s in seq_len(G)) {
    for (r in seq_len(G)) {
      if (s == r) next
      records[[k]] <- bm_2023_bilateral_source(io, s, r)
      k <- k + 1L
    }
  }

  out <- do.call(rbind, records)
  rownames(out) <- NULL
  out
}

#' BM_2023 pure bilateral (/sr) decomposition for all pairs
#'
#' Applies \code{bm_2023_bilateral_pure()} to all ordered pairs (s,r), s ≠ r.
#'
#' @param io bm_io object
#'
#' @return data.frame with one row per exporter-importer pair:
#'   exporter, importer, DVA_star_sr, DDC_star_sr, FVA_star_sr, FDC_star_sr, EX_sr
#' @export
bm_2023_bilateral_pure_all <- function(io) {
  stopifnot(inherits(io, "bm_io"))
  G <- io$G

  records <- list()
  k <- 1L

  for (s in seq_len(G)) {
    for (r in seq_len(G)) {
      if (s == r) next
      records[[k]] <- bm_2023_bilateral_pure(io, s, r)
      k <- k + 1L
    }
  }

  out <- do.call(rbind, records)
  rownames(out) <- NULL
  out
}

#' BM_2023 sink-based bilateral decomposition for all pairs
#'
#' Applies \code{bm_2023_bilateral_sink()} to all ordered pairs (s,r), s ≠ r.
#'
#' @param io bm_io object
#'
#' @return data.frame with one row per exporter-importer pair:
#'   exporter, importer, DVAsink_sr, DDCsink_sr, FVAsink_sr, FDCsink_sr, EX_sr
#' @export
bm_2023_bilateral_sink_all <- function(io) {
  stopifnot(inherits(io, "bm_io"))
  G <- io$G

  records <- list()
  k <- 1L

  for (s in seq_len(G)) {
    for (r in seq_len(G)) {
      if (s == r) next
      records[[k]] <- bm_2023_bilateral_sink(io, s, r)
      k <- k + 1L
    }
  }

  out <- do.call(rbind, records)
  rownames(out) <- NULL
  out
}

