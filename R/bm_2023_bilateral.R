# R/bm_2023_bilateral.R

#' BM_2023 source-based bilateral decomposition of exports from s to r
#'
#' @param io A \code{bm_io} object.
#' @param s Exporter country (name or index).
#' @param r Importer country (name or index).
#' @return A data frame with the source-based value-added decomposition.
#' @export
bm_2023_bilateral_source <- function(io, s, r) {
  stopifnot(inherits(io, "bm_io"))
  s <- bm_country_id(io, s)
  r <- bm_country_id(io, r)
  if (s == r) stop("s and r must be different")

  G    <- io$G
  N    <- io$N
  vvec <- io$v
  A    <- io$A
  B    <- io$B

  idx_s <- bm_idx_country(io, s)
  e_sr <- bm_get_e_sr(io, s, r)
  v_s  <- vvec[idx_s]

  Bs  <- bm_Bslash_s(io, s)
  Bss <- bm_block(io, Bs, s, s)

  sum_AB_e <- rep(0, N)
  for (j in seq_len(G)) {
    if (j == s) next
    idx_j <- bm_idx_country(io, j)
    A_sj  <- bm_block(io, A, s, j)
    B_js  <- bm_block(io, B, j, s)
    sum_AB_e <- sum_AB_e + A_sj %*% (B_js %*% e_sr)
  }

  vec_DVA   <- Bss %*% e_sr
  DVAsource <- sum(v_s * as.numeric(vec_DVA))

  vec_DDC   <- Bss %*% sum_AB_e
  DDCsource <- sum(v_s * as.numeric(vec_DDC))

  FVAsource <- 0
  FDCsource <- 0
  for (t in seq_len(G)) {
    if (t == s) next
    idx_t <- bm_idx_country(io, t)
    v_t   <- vvec[idx_t]

    Bts_slash <- bm_block(io, Bs, t, s)

    vec_FVA_t <- Bts_slash %*% e_sr
    FVAsource <- FVAsource + sum(v_t * as.numeric(vec_FVA_t))

    vec_FDC_t <- Bts_slash %*% sum_AB_e
    FDCsource <- FDCsource + sum(v_t * as.numeric(vec_FDC_t))
  }

  EX_sr <- sum(e_sr)

  data.frame(
    exporter      = io$countries[s],
    importer      = io$countries[r],
    DVAsource_sr  = DVAsource,
    DDCsource_sr  = DDCsource,
    FVAsource_sr  = FVAsource,
    FDCsource_sr  = FDCsource,
    EX_sr         = EX_sr
  )
}

#' BM_2023 pure bilateral decomposition of exports from s to r
#'
#' @param io A \code{bm_io} object.
#' @param s Exporter country (name or index).
#' @param r Importer country (name or index).
#' @return A data frame with the pure bilateral value-added decomposition.
#' @export
bm_2023_bilateral_pure <- function(io, s, r) {
  stopifnot(inherits(io, "bm_io"))
  s <- bm_country_id(io, s)
  r <- bm_country_id(io, r)
  if (s == r) stop("s and r must be different")

  G    <- io$G
  GN   <- io$GN
  vvec <- io$v
  A    <- io$A
  B    <- io$B

  idx_s <- bm_idx_country(io, s)
  v_s      <- vvec[idx_s]
  e_sr_vec <- bm_get_e_sr(io, s, r)

  A_sr <- bm_block(io, A, s, r)
  B_rs <- bm_block(io, B, r, s)

  A_mod  <- bm_Aslash_sr(io, s, r)
  I_GN   <- Matrix::Diagonal(n = GN)
  Bsr    <- solve(I_GN - A_mod)

  Bsr_ss <- bm_block(io, Bsr, s, s)

  vec_DVA  <- Bsr_ss %*% e_sr_vec
  DVA_star <- sum(v_s * as.numeric(vec_DVA))

  inter_path <- A_sr %*% (B_rs %*% e_sr_vec)

  vec_DDC  <- Bsr_ss %*% inter_path
  DDC_star <- sum(v_s * as.numeric(vec_DDC))

  FVA_star <- 0
  FDC_star <- 0
  for (t in seq_len(G)) {
    if (t == s) next
    idx_t <- bm_idx_country(io, t)
    v_t   <- vvec[idx_t]

    Bsr_ts <- bm_block(io, Bsr, t, s)

    vec_FVA_t <- Bsr_ts %*% e_sr_vec
    FVA_star  <- FVA_star + sum(v_t * as.numeric(vec_FVA_t))

    vec_FDC_t <- Bsr_ts %*% inter_path
    FDC_star  <- FDC_star + sum(v_t * as.numeric(vec_FDC_t))
  }

  EX_sr <- sum(e_sr_vec)

  data.frame(
    exporter    = io$countries[s],
    importer    = io$countries[r],
    DVA_star_sr = DVA_star,
    DDC_star_sr = DDC_star,
    FVA_star_sr = FVA_star,
    FDC_star_sr = FDC_star,
    EX_sr       = EX_sr
  )
}
