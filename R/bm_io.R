# R/bm_io.R

#' Build a BM-style IO object
#'
#' @param Z  GN x GN intermediate-use matrix
#' @param Y  GN x G  final demand by destination country
#' @param VA length-GN vector of value added
#' @param X  length-GN vector of gross output
#' @param countries character vector of country codes/names (length G)
#' @param sectors   character vector of sector codes/names (length N)
#'
#' @return A list of IO objects and metadata used by bmGVC
#' @export
bm_build_io <- function(Z, Y, VA, X, countries, sectors) {
  stopifnot(is.matrix(Z), is.matrix(Y))
  G  <- length(countries)
  N  <- length(sectors)
  GN <- G * N
  stopifnot(nrow(Z) == GN, ncol(Z) == GN)
  stopifnot(nrow(Y) == GN, ncol(Y) == G)
  stopifnot(length(VA) == GN, length(X) == GN)

  # avoid division by zero
  X_safe <- X
  X_safe[X_safe == 0] <- 1e-10

  X_hat <- Matrix::Diagonal(x = X_safe)

  # technical coefficients
  A <- Z %*% solve(X_hat)
  I_GN <- Matrix::Diagonal(n = GN)
  B <- solve(I_GN - A)

  # value-added coefficients (1 x GN row, numeric vector)
  v_row <- matrix(VA, nrow = 1) %*% solve(X_hat)
  vvec  <- as.numeric(v_row)

  # simple index helpers inside object (by position)
  idx_country <- function(g) {
    ((g - 1L) * N + 1L):(g * N)
  }

  structure(
    list(
      Z = Z,
      Y = Y,
      VA = VA,
      X  = X_safe,
      A  = A,
      B  = B,
      v  = vvec,
      G  = G,
      N  = N,
      GN = GN,
      countries      = countries,
      sectors        = sectors,
      country_codes  = stats::setNames(seq_len(G), countries),
      sector_codes   = stats::setNames(seq_len(N), sectors),
      .idx_country   = idx_country
    ),
    class = "bm_io"
  )
}
