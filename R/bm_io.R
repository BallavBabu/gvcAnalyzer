# bm_io.R
# Core IO object and helpers for bmGVC

#' Build bm_io object from MRIO core blocks
#'
#' @param Z Intermediate input matrix (GN x GN).
#'   Rows are supplying country–sector combinations.
#'   Columns are using country–sector combinations.
#' @param Y Final demand matrix (GN x G).
#'   Rows are producing country–sector combinations.
#'   Columns are destination countries.
#' @param VA Value added vector (length GN).
#'   Same row order as Z.
#' @param X Gross output vector (length GN).
#'   Same row order as Z.
#' @param countries Character vector of country names (length G).
#' @param sectors Character vector of sector names (length N).
#'
#' @return An object of class "bm_io" containing IO matrices and helpers.
#' @export
bm_build_io <- function(Z, Y, VA, X, countries, sectors) {
  stopifnot(is.matrix(Z), is.matrix(Y))
  stopifnot(is.numeric(VA), is.numeric(X))

  G  <- length(countries)
  N  <- length(sectors)
  GN <- G * N

  # Basic dimension checks
  stopifnot(nrow(Z) == GN, ncol(Z) == GN)
  stopifnot(length(VA) == GN, length(X) == GN)
  stopifnot(nrow(Y) == GN, ncol(Y) == G)

  # Avoid division by zero in coefficient construction
  X[X == 0] <- 1e-10
  X_hat <- Matrix::Diagonal(x = X)

  # Technical coefficients and global Leontief inverse
  A    <- Z %*% solve(X_hat)
  I_GN <- Matrix::Diagonal(n = GN)
  B    <- solve(I_GN - A)

  # Value added coefficients v (1 x GN), stored as numeric vector
  v_row <- matrix(VA, nrow = 1) %*% solve(X_hat)
  vvec  <- as.numeric(v_row)

  # Domestic Leontief inverse for each country g
  L_list <- vector("list", G)
  for (g in seq_len(G)) {
    idx_g <- bm_idx_country(g, N)
    A_gg  <- A[idx_g, idx_g, drop = FALSE]
    L_list[[g]] <- solve(Matrix::Diagonal(N) - A_gg)
  }

  structure(
    list(
      Z         = Z,
      Y         = Y,
      VA        = VA,
      X         = X,
      countries = countries,
      sectors   = sectors,
      G         = G,
      N         = N,
      GN        = GN,
      A         = A,
      B         = B,
      v         = vvec,
      L_list    = L_list
    ),
    class = "bm_io"
  )
}

# Internal helper: indices of country g (1..G) in stacked GN vector
bm_idx_country <- function(g, N) {
  ((g - 1L) * N + 1L):(g * N)
}

# Internal helper: N x N block (i,j) from GN x GN matrix M
bm_block <- function(M, i, j, N) {
  idx_i <- bm_idx_country(i, N)
  idx_j <- bm_idx_country(j, N)
  M[idx_i, idx_j, drop = FALSE]
}

#' Bilateral exports e_sr by sector
#'
#' Total exports from country s to country r by sector, including
#' both intermediate and final goods.
#'
#' @param io bm_io object.
#' @param s Exporter index (1..G).
#' @param r Importer index (1..G).
#'
#' @return Numeric vector of length N with exports from s to r by sector.
#' @export
bm_get_e_sr <- function(io, s, r) {
  G  <- io$G
  N  <- io$N
  Z  <- io$Z
  Y  <- io$Y

  stopifnot(s %in% seq_len(G), r %in% seq_len(G))

  idx_s <- bm_idx_country(s, N)
  idx_r <- bm_idx_country(r, N)

  int_sr <- rowSums(Z[idx_s, idx_r, drop = FALSE])  # intermediates s -> r
  fin_sr <- Y[idx_s, r]                             # final goods s -> r

  int_sr + fin_sr
}

#' Country-total exports e_s* by sector
#'
#' Total exports from country s to all foreign destinations,
#' decomposed by sector.
#'
#' @param io bm_io object.
#' @param s Exporter index (1..G).
#'
#' @return Numeric vector of length N with exports from s to the world.
#' @export
bm_get_e_star <- function(io, s) {
  G  <- io$G
  N  <- io$N
  GN <- io$GN
  Z  <- io$Z
  Y  <- io$Y

  stopifnot(s %in% seq_len(G))

  idx_s <- bm_idx_country(s, N)

  # Intermediate exports to all foreign country–sectors
  int_exp <- rowSums(Z[idx_s, setdiff(seq_len(GN), idx_s), drop = FALSE])

  # Final exports to all foreign destinations
  fin_exp <- rowSums(Y[idx_s, setdiff(seq_len(G), s), drop = FALSE])

  int_exp + fin_exp
}
