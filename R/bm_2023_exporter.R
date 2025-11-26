# R/bm_2023_exporter.R

#' BM_2023 exporter-perspective decomposition of total exports of s
#'
#' @param io A \code{bm_io} object.
#' @param s Exporter country (name or index).
#' @return A data frame with the exporter-total decomposition.
#' @export
bm_2023_exporter_total <- function(io, s) {
  stopifnot(inherits(io, "bm_io"))
  s <- bm_country_id(io, s)

  G  <- io$G
  N  <- io$N
  GN <- io$GN

  idx_s <- bm_idx_country(io, s)

  vvec <- io$v
  v_s  <- vvec[idx_s]

  A    <- io$A
  B    <- io$B

  e_s  <- bm_get_e_star(io, s)
  Bs   <- bm_Bslash_s(io, s)
  Bss  <- bm_block(io, Bs, s, s)

  DVA_s <- sum(v_s * as.numeric(Bss %*% e_s))

  sum_AB_e <- rep(0, N)
  for (j in seq_len(G)) {
    if (j == s) next
    A_sj  <- bm_block(io, A, s, j)
    B_js  <- bm_block(io, B, j, s)
    sum_AB_e <- sum_AB_e + A_sj %*% (B_js %*% e_s)
  }

  DDC_s <- sum(v_s * as.numeric(Bss %*% sum_AB_e))

  FVA_s <- 0
  FDC_s <- 0
  for (t in seq_len(G)) {
    if (t == s) next
    idx_t <- bm_idx_country(io, t)
    v_t   <- vvec[idx_t]

    Bts <- bm_block(io, Bs, t, s)

    vec_FVA_t <- Bts %*% e_s
    FVA_s <- FVA_s + sum(v_t * as.numeric(vec_FVA_t))

    vec_FDC_t <- Bts %*% sum_AB_e
    FDC_s <- FDC_s + sum(v_t * as.numeric(vec_FDC_t))
  }

  EX_s <- sum(e_s)

  data.frame(
    country = io$countries[s],
    DVA_s   = DVA_s,
    DDC_s   = DDC_s,
    FVA_s   = FVA_s,
    FDC_s   = FDC_s,
    EX_s    = EX_s
  )
}
