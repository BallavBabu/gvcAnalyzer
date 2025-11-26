# R/bm_io.R

#' Build a bm_io object from IO table blocks
#'
#' @param Z Intermediate demand matrix (GN x GN).
#' @param Y Final demand matrix. Can be (GN x G) OR (GN x (G * FD_categories)).
#' @param VA Value added. Can be a vector (length GN) or matrix (Rows x GN).
#' @param X Output vector (length GN).
#' @param countries Character vector of country names/codes (length G).
#' @param sectors Character vector of sector names/codes (length N).
#'
#' @return An object of class \code{"bm_io"}.
#' @export
#' @importFrom Matrix Diagonal
bm_build_io <- function(Z, Y, VA, X, countries, sectors) {
  Z <- as.matrix(Z)
  Y <- as.matrix(Y)
  X <- as.numeric(X)

  G <- length(countries)
  N <- length(sectors)
  GN <- G * N

  if (nrow(Z) != GN || ncol(Z) != GN) stop("Z dimensions mismatch.")

  if (is.matrix(VA)) {
    if (ncol(VA) != GN) stop("VA matrix columns must equal GN.")
    VA <- colSums(VA)
  }
  VA <- as.numeric(VA)

  if (ncol(Y) == G) {
    Y_agg <- Y
  } else if (ncol(Y) %% G == 0) {
    n_fd <- ncol(Y) / G
    Y_agg <- matrix(0, nrow = GN, ncol = G)
    for (g in 1:G) {
      idx_start <- (g - 1) * n_fd + 1
      idx_end   <- g * n_fd
      Y_agg[, g] <- rowSums(Y[, idx_start:idx_end, drop = FALSE])
    }
    colnames(Y_agg) <- countries
  } else {
    stop("Y columns must be a multiple of G.")
  }

  X_adj <- X
  X_adj[X_adj == 0] <- 1e-10
  X_hat_inv <- Matrix::Diagonal(x = 1/X_adj)

  A <- Z %*% X_hat_inv
  I_GN <- Matrix::Diagonal(n = GN)
  B <- solve(I_GN - A)
  v <- as.numeric(matrix(VA, nrow = 1) %*% X_hat_inv)

  country_codes <- as.list(seq_along(countries))
  names(country_codes) <- countries

  structure(list(
    Z = Z, Y = Y_agg, VA = VA, X = X,
    countries = countries, sectors = sectors,
    country_codes = country_codes,
    G = G, N = N, GN = GN,
    A = A, B = B, v = v
  ), class = "bm_io")
}

#' Exports from s to r (sectoral)
#'
#' @param io bm_io object.
#' @param exporter Exporter country (name or index).
#' @param importer Importer country (name or index).
#' @return Numeric vector of exports.
#' @export
bm_get_e_sr <- function(io, exporter, importer) {
  s_id <- bm_country_id(io, exporter)
  r_id <- bm_country_id(io, importer)
  idx_s <- bm_idx_country(io, s_id)
  idx_r <- bm_idx_country(io, r_id)
  rowSums(io$Z[idx_s, idx_r, drop = FALSE]) + io$Y[idx_s, r_id]
}

#' Total exports of s to all foreign destinations
#'
#' @param io bm_io object.
#' @param exporter Exporter country (name or index).
#' @return Numeric vector of total exports.
#' @export
bm_get_e_star <- function(io, exporter) {
  s_id <- bm_country_id(io, exporter)
  idx_s <- bm_idx_country(io, s_id)
  idx_foreign_Z <- setdiff(seq_len(io$GN), idx_s)
  idx_foreign_Y <- setdiff(seq_len(io$G), s_id)
  rowSums(io$Z[idx_s, idx_foreign_Z, drop = FALSE]) +
    rowSums(io$Y[idx_s, idx_foreign_Y, drop = FALSE])
}
