# R/bm_2023_bilateral.R
#
# BM_2023 bilateral decompositions:
#   1) Source-based (/s) exporter perspective
#   2) Pure bilateral (/sr) perspective

# -------------------------------------------------------------------
# Internal helpers
# -------------------------------------------------------------------

# Country index range  (delegates to bm_io internal)
bm_idx_country <- function(io, g) {
  io$.idx_country(g)
}

# Block (j,k) of a slash-s Leontief Bs = B^/s
bm_block_slash_s <- function(io, Bs, j, k) {
  Bs[bm_idx_country(io, j), bm_idx_country(io, k), drop = FALSE]
}

# -------------------------------------------------------------------
# 1) Source-based bilateral BM_2023 decomposition (/s)
# -------------------------------------------------------------------

#' BM_2023 source-based bilateral decomposition of exports from s to r
#'
#' Decomposes gross exports e_{sr} (from exporter s to importer r) as
#'   u_N × e_{sr} = DVAsource_sr + DDCsource_sr + FVAsource_sr + FDCsource_sr
#' using the exporter-based slash matrix A^/s (and B^/s).
#'
#' @param io bm_io object as returned by \code{bm_build_io()}
#' @param s  exporter (country index or country code, e.g. 1 or "China")
#' @param r  importer (country index or country code, e.g. 2 or "India")
#'
#' @return A data frame with one row for the pair (s,r):
#'   \itemize{
#'     \item \code{DVAsource_sr} Domestic value added from s in e_{sr}
#'           that never returns to s.
#'     \item \code{DDCsource_sr} Double-counting term from s due to
#'           multi-border domestic returns in e_{sr}.
#'     \item \code{FVAsource_sr} Foreign value added in e_{sr}.
#'     \item \code{FDCsource_sr} Foreign double-counting in e_{sr}.
#'     \item \code{EX_sr}        Gross exports e_{sr} (scalar).
#'   }
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

  # Bilateral exports e_sr (N-vector) and VA coeffs for s
  e_sr <- bm_e_sr(io, s, r)
  v_s  <- vvec[idx_s]

  # B^/s and its (s,s) block
  Bs  <- bm_Bslash_s(io, s)              # GN x GN, slash-s Leontief
  Bss <- bm_block_slash_s(io, Bs, s, s)  # N x N

  # Σ_{j≠s}(A_sj B_js e_sr)
  sum_AB_e <- rep(0, N)
  for (j in seq_len(G)) {
    if (j == s) next
    idx_j <- bm_idx_country(io, j)
    A_sj  <- A[idx_s, idx_j, drop = FALSE]
    B_js  <- B[idx_j, idx_s, drop = FALSE]
    sum_AB_e <- sum_AB_e + A_sj %*% (B_js %*% e_sr)
  }

  # DVAsource_sr = v_s × B^/s_ss × e_sr
  vec_DVA   <- Bss %*% e_sr
  DVAsource <- sum(v_s * as.numeric(vec_DVA))

  # DDCsource_sr = v_s × B^/s_ss × Σ_{j≠s}(A_sj B_js e_sr)
  vec_DDC   <- Bss %*% sum_AB_e
  DDCsource <- sum(v_s * as.numeric(vec_DDC))

  # FVAsource_sr = Σ_{t≠s} v_t × B^/s_ts × e_sr
  # FDCsource_sr = Σ_{t≠s} v_t × B^/s_ts × Σ_{j≠s}(A_sj B_js e_sr)
  FVAsource <- 0
  FDCsource <- 0
  for (t in seq_len(G)) {
    if (t == s) next
    idx_t <- bm_idx_country(io, t)
    v_t   <- vvec[idx_t]

    Bts_slash <- bm_block_slash_s(io, Bs, t, s)  # N x N

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

# -------------------------------------------------------------------
# 2) Pure bilateral /sr BM_2023 decomposition
# -------------------------------------------------------------------

# Internal: A^/sr = A with the bilateral intermediate block A_sr set to zero
bm_Aslash_sr <- function(io, s, r) {
  A_mod  <- io$A
  rows_s <- bm_idx_country(io, s)
  cols_r <- bm_idx_country(io, r)
  A_mod[rows_s, cols_r] <- 0
  A_mod
}

#' BM_2023 pure bilateral decomposition of exports from s to r
#'
#' Implements the /sr-based decomposition:
#'   u_N × e_{sr} = DVA_star_sr + DDC_star_sr + FVA_star_sr + FDC_star_sr
#' where A^/sr zeros only the bilateral intermediate block A_{sr}.
#'
#' @param io bm_io object as returned by \code{bm_build_io()}
#' @param s  exporter (country index or code)
#' @param r  importer (country index or code)
#'
#' @return A data frame with one row for the pair (s,r):
#'   \itemize{
#'     \item \code{DVA_star_sr} Domestic VA of s in bilateral exports e_{sr},
#'           when only the A_{sr} block is cut.
#'     \item \code{DDC_star_sr} Double-counting term for s in e_{sr}.
#'     \item \code{FVA_star_sr} Foreign VA in e_{sr}.
#'     \item \code{FDC_star_sr} Foreign double-counting in e_{sr}.
#'     \item \code{EX_sr}       Gross exports e_{sr} (scalar).
#'   }
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
  idx_r <- bm_idx_country(io, r)

  v_s      <- vvec[idx_s]
  e_sr_vec <- bm_e_sr(io, s, r)                 # N-vector

  # Original blocks
  A_sr <- A[idx_s, idx_r, drop = FALSE]        # N x N
  B_rs <- B[idx_r, idx_s, drop = FALSE]        # N x N

  # Bilateral slash Leontief B^/sr = (I - A^/sr)^(-1)
  A_mod  <- bm_Aslash_sr(io, s, r)
  I_GN   <- Matrix::Diagonal(n = GN)
  Bsr    <- solve(I_GN - A_mod)

  # (s,s) block of B^/sr
  Bsr_ss <- Bsr[idx_s, idx_s, drop = FALSE]

  # DVA*_sr = v_s × B^/sr_ss × e_sr
  vec_DVA  <- Bsr_ss %*% e_sr_vec
  DVA_star <- sum(v_s * as.numeric(vec_DVA))

  # Intermediate path via A_sr × B_rs × e_sr
  inter_path <- A_sr %*% (B_rs %*% e_sr_vec)

  # DDC*_sr = v_s × B^/sr_ss × A_sr × B_rs × e_sr
  vec_DDC  <- Bsr_ss %*% inter_path
  DDC_star <- sum(v_s * as.numeric(vec_DDC))

  # Foreign terms: FVA*_sr and FDC*_sr
  FVA_star <- 0
  FDC_star <- 0
  for (t in seq_len(G)) {
    if (t == s) next
    idx_t <- bm_idx_country(io, t)
    v_t   <- vvec[idx_t]

    Bsr_ts <- Bsr[idx_t, idx_s, drop = FALSE]

    # FVA* contribution from origin t
    vec_FVA_t <- Bsr_ts %*% e_sr_vec
    FVA_star  <- FVA_star + sum(v_t * as.numeric(vec_FVA_t))

    # FDC* contribution from origin t
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
