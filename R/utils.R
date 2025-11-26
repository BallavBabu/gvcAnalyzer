# R/utils.R

#' @importFrom stats aggregate
#' @import methods
NULL

#' Internal: Resolve country to index
#' @noRd
bm_country_id <- function(io, g) {
  if (is.numeric(g)) return(as.integer(g))
  if (is.character(g)) {
    idx <- io$country_codes[[g]]
    if (is.null(idx)) stop('Country not found: ', g)
    return(idx)
  }
  stop('Invalid country identifier')
}

#' Internal: Get indices for country g in GN matrix
#' @noRd
bm_idx_country <- function(io, g) {
  id <- bm_country_id(io, g)
  ((id - 1L) * io$N + 1L):(id * io$N)
}

#' Internal: Get block (s, r) from matrix M
#' @noRd
bm_block <- function(io, M, s, r) {
  idx_s <- bm_idx_country(io, s)
  idx_r <- bm_idx_country(io, r)
  M[idx_s, idx_r, drop = FALSE]
}

#' Internal: Precompute domestic Leontief inverses
#' @noRd
bm_L_list <- function(io) {
  L_list <- vector('list', io$G)
  I_N <- Matrix::Diagonal(n = io$N)
  for (g in 1:io$G) {
    idx <- bm_idx_country(io, g)
    A_gg <- io$A[idx, idx, drop = FALSE]
    L_list[[g]] <- solve(I_N - A_gg)
  }
  L_list
}

# --- BM 2023 Legacy Helpers (Slash Matrices) ---

#' Internal: A^/s
#' @noRd
bm_Aslash_s <- function(io, s) {
  A <- io$A
  G <- io$G
  idx_s <- bm_idx_country(io, s)
  for (j in seq_len(G)) {
    if (j == bm_country_id(io, s)) next
    idx_j <- bm_idx_country(io, j)
    A[idx_s, idx_j] <- 0
  }
  A
}

#' Internal: B^/s
#' @noRd
bm_Bslash_s <- function(io, s) {
  GN <- io$GN
  As <- bm_Aslash_s(io, s)
  I_GN <- Matrix::Diagonal(n = GN)
  solve(I_GN - As)
}

#' Internal: A^/sr
#' @noRd
bm_Aslash_sr <- function(io, s, r) {
  A_mod  <- io$A
  rows_s <- bm_idx_country(io, s)
  cols_r <- bm_idx_country(io, r)
  A_mod[rows_s, cols_r] <- 0
  A_mod
}
