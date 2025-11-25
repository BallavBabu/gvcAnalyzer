# R/bm_io.R
#
# Core IO constructor and helpers for bmGVC
# - bm_build_io(): builds bm_io object with A, B, v, etc.
# - bm_get_e_sr(): exports from s to r by sector (intermediate + final)
# - bm_get_e_star(): exports from s to all foreign partners by sector

#' Build a bm_io object from IO table blocks
#'
#' @param Z Intermediate demand matrix (GN x GN), rows = suppliers, cols = users.
#' @param Y Final demand matrix (GN x G), columns = destination countries.
#' @param VA Value added vector (length GN), in same industry order as Z rows.
#' @param X Output vector (length GN), in same industry order as Z rows.
#' @param countries Character vector of country names, length G.
#' @param sectors Character vector of sector names, length N.
#'
#' @return An object of class \code{"bm_io"} containing
#'   \itemize{
#'     \item \code{Z, Y, VA, X}
#'     \item \code{countries, sectors}
#'     \item \code{G, N, GN}
#'     \item \code{A}: technical coefficients (GN x GN)
#'     \item \code{B}: Leontief inverse (GN x GN)
#'     \item \code{v}: value-added coefficients (length GN)
#'   }
#' @export
bm_build_io <- function(Z, Y, VA, X, countries, sectors) {
  # Basic dimensions
  Z <- as.matrix(Z)
  Y <- as.matrix(Y)
  VA <- as.numeric(VA)
  X  <- as.numeric(X)

  GN <- nrow(Z)
  if (ncol(Z) != GN) {
    stop("Z must be square GN x GN.")
  }
  if (length(VA) != GN || length(X) != GN) {
    stop("VA and X must have length equal to nrow(Z).")
  }

  G <- length(countries)
  N <- length(sectors)
  if (GN != G * N) {
    stop("GN != G*N. Check countries/sectors vs Z dimension.")
  }

  if (!all(dim(Y) == c(GN, G))) {
    stop("Y must be GN x G (rows = industries, cols = destination countries).")
  }

  # Avoid division by zero
  X_adj <- X
  X_adj[X_adj == 0] <- 1e-10

  X_hat <- Matrix::Diagonal(x = X_adj)

  # Technical coefficients and Leontief inverse
  A <- Z %*% solve(X_hat)
  I_GN <- Matrix::Diagonal(n = GN)
  B <- solve(I_GN - A)

  # Value-added coefficients
  v_row <- matrix(VA, nrow = 1) %*% solve(X_hat)  # 1 x GN
  vvec  <- as.numeric(v_row)

  out <- list(
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
    v         = vvec
  )
  class(out) <- "bm_io"
  out
}

# ------------------------------------------------------------
# Internal helpers for bm_io
# ------------------------------------------------------------

bm_country_id <- function(io, g) {
  # g can be numeric index or country name
  if (is.character(g)) {
    idx <- match(g, io$countries)
    if (is.na(idx)) stop("Unknown country name: ", g)
    return(idx)
  }
  if (!is.numeric(g) || length(g) != 1L) {
    stop("Invalid country identifier.")
  }
  if (g < 1L || g > io$G) {
    stop("Country index out of range.")
  }
  as.integer(g)
}

bm_idx_country <- function(io, g) {
  g_id <- bm_country_id(io, g)
  ((g_id - 1L) * io$N + 1L):(g_id * io$N)
}

bm_block_io <- function(io, M, g_row, g_col) {
  # extract N x N block of matrix M corresponding to g_row, g_col
  idx_r <- bm_idx_country(io, g_row)
  idx_c <- bm_idx_country(io, g_col)
  M[idx_r, idx_c, drop = FALSE]
}

# ------------------------------------------------------------
# Export helpers: e_sr and e_s*
# ------------------------------------------------------------

#' Exports from s to r by sector (intermediate + final)
#'
#' Computes a length-N vector of exports from exporter \eqn{s} to importer \eqn{r},
#' summing intermediate and final exports:
#' \deqn{ e_{sr} = Z_{s,r} \mathbf{u}_N + Y_{s,r} }
#'
#' \code{exporter} and \code{importer} can be numeric indices or country names.
#'
#' @param io bm_io object
#' @param exporter exporter country (index or name)
#' @param importer importer country (index or name)
#'
#' @return Numeric vector of length N with sectoral exports s -> r.
#' @export
bm_get_e_sr <- function(io, exporter, importer) {
  stopifnot(inherits(io, "bm_io"))

  s_id <- bm_country_id(io, exporter)
  r_id <- bm_country_id(io, importer)

  rows_s <- bm_idx_country(io, s_id)
  cols_r <- bm_idx_country(io, r_id)

  # intermediate exports: Z_sr * u
  int_sr <- rowSums(io$Z[rows_s, cols_r, drop = FALSE])

  # final exports: Y_sr (column r of Y)
  fin_sr <- io$Y[rows_s, r_id]

  int_sr + fin_sr
}

#' Exports from s to all foreign partners by sector (e_s*)
#'
#' Computes a length-N vector of total exports from exporter \eqn{s} to all
#' foreign partners (intermediate + final), excluding domestic uses:
#' \deqn{ e_s^* = \sum_{r \neq s} e_{sr} }
#'
#' \code{exporter} can be numeric index or country name.
#'
#' @param io bm_io object
#' @param exporter exporter country (index or name)
#'
#' @return Numeric vector of length N with sectoral exports of s to all foreign partners.
#' @export
bm_get_e_star <- function(io, exporter) {
  stopifnot(inherits(io, "bm_io"))

  s_id <- bm_country_id(io, exporter)

  rows_s <- bm_idx_country(io, s_id)
  # intermediate exports: sum over all foreign columns
  foreign_cols <- setdiff(seq_len(io$GN), rows_s)
  int_exp <- rowSums(io$Z[rows_s, foreign_cols, drop = FALSE])

  # final exports: sum over all foreign final-demand columns
  foreign_fd <- setdiff(seq_len(io$G), s_id)
  fin_exp <- rowSums(io$Y[rows_s, foreign_fd, drop = FALSE])

  int_exp + fin_exp
}
