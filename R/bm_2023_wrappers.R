# R/bm_2023_wrappers.R

#' BM_2023 exporter totals for all countries
#' @param io A \code{bm_io} object.
#' @return Data frame of exporter totals for all countries.
#' @export
bm_2023_exporter_total_all <- function(io) {
  stopifnot(inherits(io, "bm_io"))
  do.call(
    rbind,
    lapply(seq_len(io$G), function(s) bm_2023_exporter_total(io, s))
  )
}

#' BM_2023 source-based bilateral decomposition for all pairs
#' @param io A \code{bm_io} object.
#' @return Data frame of source-based decomposition for all pairs.
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
  do.call(rbind, records)
}

#' BM_2023 pure bilateral (/sr) decomposition for all pairs
#' @param io A \code{bm_io} object.
#' @return Data frame of pure bilateral decomposition for all pairs.
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
  do.call(rbind, records)
}

#' BM_2023 sink-based bilateral decomposition for all pairs
#' @param io A \code{bm_io} object.
#' @return Data frame of sink-based decomposition for all pairs.
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
  do.call(rbind, records)
}
