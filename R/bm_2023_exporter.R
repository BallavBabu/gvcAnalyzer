# R/bm_2023_exporter.R

# Internal helper: index range for country g
bm_idx_country <- function(io, g) {
  io$.idx_country(g)
}

# Internal helper: block (i,j) of a GN x GN matrix
bm_block <- function(M, idx_i, idx_j) {
  M[idx_i, idx_j, drop = FALSE]
}

# Internal: total export vector e_s* (N x 1) for country s
bm_e_star <- function(io, s) {
  G  <- io$G
  N  <- io$N
  GN <- io$GN
  idx_s <- bm_idx_country(io, s)

  Z <- io$Z
  Y <- io$Y

  int_exp <- rowSums(Z[idx_s, setdiff(seq_len(GN), idx_s), drop = FALSE])
  fin_exp <- rowSums(Y[idx_s, setdiff(seq_len(G), s),   drop = FALSE])
  as.numeric(int_exp + fin_exp)
}

# Internal: bilateral export vector e_sr (N x 1) from s to r
bm_e_sr <- function(io, s, r) {
  idx_s <- bm_idx_country(io, s)
  idx_r <- bm_idx_country(io, r)

  Z <- io$Z
  Y <- io$Y

  int_sr <- rowSums(Z[idx_s, idx_r, drop = FALSE])
  fin_sr <- Y[idx_s, r]
  as.numeric(int_sr + fin_sr)
}

# Internal: build A^/s (zero all intermediate exports from s to j != s)
bm_Aslash_s <- function(io, s) {
  A <- io$A
  G <- io$G
  N <- io$N

  idx_s <- bm_idx_country(io, s)

  for (j in seq_len(G)) {
    if (j == s) next
    idx_j <- bm_idx_country(io, j)
    A[idx_s, idx_j] <- 0
  }
  A
}

# Internal: B^/s = (I - A^/s)^(-1)
bm_Bslash_s <- function(io, s) {
  GN <- io$GN
  As <- bm_Aslash_s(io, s)
  I_GN <- Matrix::Diagonal(n = GN)
  solve(I_GN - As)
}

#' BM_2023 exporter-perspective decomposition of total exports of s
#'
#' Implements Eq. (4) at the country level:
#' u_N × e_s* = DVA_s* + DDC_s* + FVA_s* + FDC_s*
#'
#' @param io bm_io object
#' @param s  exporting country (index or code, e.g. 1 or "China")
#'
#' @return data.frame with one row for country s
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

  e_s  <- bm_e_star(io, s)          # N x 1
  Bs   <- bm_Bslash_s(io, s)        # GN x GN
  Bss  <- bm_block(Bs, idx_s, idx_s)

  # ----- skeleton for terms, to be filled with exact BM_2023 formulas -----
  # DVAsource_s* = v_s × B^/s_ss × e_s*
  DVA_s <- sum(v_s * as.numeric(Bss %*% e_s))

  # Σ_{j≠s}(A_sj B_js e_s*)
  sum_AB_e <- rep(0, N)
  for (j in seq_len(G)) {
    if (j == s) next
    idx_j <- bm_idx_country(io, j)
    A_sj  <- bm_block(A, idx_s, idx_j)
    B_js  <- bm_block(B, idx_j, idx_s)
    sum_AB_e <- sum_AB_e + A_sj %*% (B_js %*% e_s)
  }

  # DDC_s* = v_s × B^/s_ss × Σ_{j≠s}(A_sj B_js e_s*)
  DDC_s <- sum(v_s * as.numeric(Bss %*% sum_AB_e))

  # FVA_s* = Σ_{t≠s} v_t × B^/s_ts × e_s*
  FVA_s <- 0
  FDC_s <- 0
  for (t in seq_len(G)) {
    if (t == s) next
    idx_t <- bm_idx_country(io, t)
    v_t   <- vvec[idx_t]

    Bts <- bm_block(Bs, idx_t, idx_s)

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
