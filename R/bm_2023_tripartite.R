# R/bm_2023_tripartite.R
#
# BM_2023 trade-based tripartite GVC decomposition:
#   For each exporter s and importer r:
#     E_sr     = gross exports from s to r
#     DAVAX_sr = domestic VA from s that crosses only one border
#     GVC_sr   = E_sr - DAVAX_sr
#     GVC_PF   = pure-forward GVC trade
#     GVC_TS   = two-sided GVC trade
#     GVC_PB   = pure-backward GVC trade
#
# All functions accept exporter/importer as index (1,2,...) or name ("China", ...).

# ------------------------------------------------------------
# Local helpers
# ------------------------------------------------------------

bm23_idx_country <- function(io, g) {
  # g can be numeric or country name
  if (is.character(g)) {
    g <- match(g, io$countries)
  }
  if (is.na(g) || g < 1L || g > io$G) {
    stop("Invalid country identifier in bm23_idx_country().")
  }
  ((g - 1L) * io$N + 1L):(g * io$N)
}

bm23_L_list <- function(io) {
  G <- io$G
  N <- io$N
  A <- io$A

  L_list <- vector("list", G)
  for (g in seq_len(G)) {
    idx_g <- bm23_idx_country(io, g)
    A_gg  <- as.matrix(A[idx_g, idx_g, drop = FALSE])
    L_list[[g]] <- solve(Matrix::Diagonal(N) - A_gg)
  }
  L_list
}

# ------------------------------------------------------------
# Core bilateral building blocks for BM_2023
# ------------------------------------------------------------

# E_sr: gross exports (intermediate + final) from s to r
bm23_E_sr <- function(io, s, r) {
  # use existing exported helper
  e_sr_vec <- bm_get_e_sr(io, s, r)
  sum(e_sr_vec)
}

# DAVAX_sr: domestic VA from s that crosses one border and
# is absorbed in r's final demand of s and r goods
bm23_DAVAX_sr <- function(io, s, r, L_list) {
  if (is.character(s)) s <- match(s, io$countries)
  if (is.character(r)) r <- match(r, io$countries)

  N   <- io$N
  A   <- io$A
  Y   <- io$Y
  v   <- io$v

  idx_s <- bm23_idx_country(io, s)
  idx_r <- bm23_idx_country(io, r)

  v_s  <- v[idx_s]                             # length N
  L_ss <- L_list[[s]]                          # N x N
  L_rr <- L_list[[r]]                          # N x N
  A_sr <- as.matrix(A[idx_s, idx_r, drop = FALSE])  # N x N
  Y_sr <- matrix(Y[idx_s, r], ncol = 1)        # N x 1
  Y_rr <- matrix(Y[idx_r, r], ncol = 1)        # N x 1

  term <- Y_sr + A_sr %*% (L_rr %*% Y_rr)      # N x 1

  as.numeric(matrix(v_s, nrow = 1) %*% (L_ss %*% term))
}

# PB_sr: pure-backward GVC trade (import content of exports s->r
# that is absorbed in r's final demand)
bm23_PB_sr <- function(io, s, r, L_list) {
  if (is.character(s)) s <- match(s, io$countries)
  if (is.character(r)) r <- match(r, io$countries)

  G   <- io$G
  N   <- io$N
  A   <- io$A
  Y   <- io$Y

  idx_s <- bm23_idx_country(io, s)
  idx_r <- bm23_idx_country(io, r)

  L_ss <- L_list[[s]]
  L_rr <- L_list[[r]]
  A_sr <- as.matrix(A[idx_s, idx_r, drop = FALSE])
  Y_sr <- matrix(Y[idx_s, r], ncol = 1)
  Y_rr <- matrix(Y[idx_r, r], ncol = 1)

  # output in s to support r's final demand of s and r goods
  term_final <- Y_sr + A_sr %*% (L_rr %*% Y_rr)  # N x 1
  q_s <- L_ss %*% term_final                     # N x 1

  u_N <- matrix(1, nrow = 1, ncol = N)

  PB_val <- 0
  for (t in seq_len(G)) {
    if (t == s) next
    idx_t <- bm23_idx_country(io, t)
    A_ts  <- as.matrix(A[idx_t, idx_s, drop = FALSE])       # inputs from t used in s
    PB_val <- PB_val + as.numeric(u_N %*% (A_ts %*% q_s))
  }
  PB_val
}

# TS_sr: two-sided GVC trade (imported inputs in s used to produce exports to r,
# which r then uses for its own exports to others)
bm23_TS_sr <- function(io, s, r, L_list) {
  if (is.character(s)) s <- match(s, io$countries)
  if (is.character(r)) r <- match(r, io$countries)

  G   <- io$G
  N   <- io$N
  A   <- io$A
  Y   <- io$Y

  idx_s <- bm23_idx_country(io, s)
  idx_r <- bm23_idx_country(io, r)

  L_ss <- L_list[[s]]
  L_rr <- L_list[[r]]
  A_sr <- as.matrix(A[idx_s, idx_r, drop = FALSE])

  # exports of r to all partners j ≠ r (by sector)
  sum_e_r_to_others <- rep(0, N)
  for (j in seq_len(G)) {
    if (j == r) next
    sum_e_r_to_others <- sum_e_r_to_others + bm_get_e_sr(io, r, j)
  }
  sum_e_r_to_others <- matrix(sum_e_r_to_others, ncol = 1)

  # demand in s arising from r's exports to others
  term_two <- A_sr %*% (L_rr %*% sum_e_r_to_others)  # N x 1
  q_s <- L_ss %*% term_two                           # N x 1

  u_N <- matrix(1, nrow = 1, ncol = N)

  TS_val <- 0
  for (t in seq_len(G)) {
    if (t == s) next
    idx_t <- bm23_idx_country(io, t)
    A_ts  <- as.matrix(A[idx_t, idx_s, drop = FALSE])
    TS_val <- TS_val + as.numeric(u_N %*% (A_ts %*% q_s))
  }
  TS_val
}

# ------------------------------------------------------------
# Exported user-level functions
# ------------------------------------------------------------

#' BM_2023 tripartite GVC trade decomposition for one pair (s,r)
#'
#' Computes BM_2023 trade-based decomposition for exporter \eqn{s} and importer \eqn{r}:
#' \itemize{
#'   \item \code{E_sr}: gross exports from \eqn{s} to \eqn{r}
#'   \item \code{DAVAX_sr}: domestic value added from \eqn{s} that crosses one border
#'   \item \code{GVC_sr}: GVC-related trade = \code{E_sr - DAVAX_sr}
#'   \item \code{GVC_PF}: pure-forward GVC trade
#'   \item \code{GVC_TS}: two-sided GVC trade
#'   \item \code{GVC_PB}: pure-backward GVC trade
#' }
#'
#' \code{exporter} and \code{importer} can be indices (1,2,...) or country names.
#'
#' @param io bm_io object
#' @param exporter exporter country (index or name)
#' @param importer importer country (index or name)
#'
#' @return data.frame with one row.
#' @export
bm_2023_tripartite_pair <- function(io, exporter, importer) {
  stopifnot(inherits(io, "bm_io"))

  # convert to indices
  if (is.character(exporter)) exporter <- match(exporter, io$countries)
  if (is.character(importer)) importer <- match(importer, io$countries)

  if (exporter == importer) {
    stop("exporter and importer must be different in bm_2023_tripartite_pair().")
  }

  L_list <- bm23_L_list(io)

  E_sr     <- bm23_E_sr(io, exporter, importer)
  DAVAX_sr <- bm23_DAVAX_sr(io, exporter, importer, L_list)
  GVC_sr   <- E_sr - DAVAX_sr

  PB_sr <- bm23_PB_sr(io, exporter, importer, L_list)
  TS_sr <- bm23_TS_sr(io, exporter, importer, L_list)
  PF_sr <- GVC_sr - PB_sr - TS_sr

  data.frame(
    exporter   = io$countries[exporter],
    importer   = io$countries[importer],
    E_sr       = E_sr,
    DAVAX_sr   = DAVAX_sr,
    GVC_sr     = GVC_sr,
    GVC_PF     = PF_sr,
    GVC_TS     = TS_sr,
    GVC_PB     = PB_sr
  )
}

#' BM_2023 tripartite GVC trade for all bilateral pairs
#'
#' Computes BM_2023 trade-based tripartite decomposition for all ordered pairs
#' \eqn{(s,r)} with \eqn{s \neq r}.
#'
#' @param io bm_io object
#'
#' @return data.frame with one row per bilateral pair.
#' @export
bm_2023_tripartite_bilateral <- function(io) {
  stopifnot(inherits(io, "bm_io"))

  G      <- io$G
  L_list <- bm23_L_list(io)

  out_list <- list()
  k <- 1L

  for (s in seq_len(G)) {
    for (r in seq_len(G)) {
      if (s == r) next

      E_sr     <- bm23_E_sr(io, s, r)
      DAVAX_sr <- bm23_DAVAX_sr(io, s, r, L_list)
      GVC_sr   <- E_sr - DAVAX_sr

      PB_sr <- bm23_PB_sr(io, s, r, L_list)
      TS_sr <- bm23_TS_sr(io, s, r, L_list)
      PF_sr <- GVC_sr - PB_sr - TS_sr

      out_list[[k]] <- data.frame(
        exporter = io$countries[s],
        importer = io$countries[r],
        E_sr     = E_sr,
        DAVAX_sr = DAVAX_sr,
        GVC_sr   = GVC_sr,
        GVC_PF   = PF_sr,
        GVC_TS   = TS_sr,
        GVC_PB   = PB_sr
      )
      k <- k + 1L
    }
  }

  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out
}

#' BM_2023 trade-based GVC components by exporter
#'
#' Aggregates BM_2023 tripartite GVC trade over all importers \eqn{r} for each exporter \eqn{s}.
#'
#' @param io bm_io object
#'
#' @return data.frame with one row per exporter and columns
#'   \code{GVC_sr, GVC_PF, GVC_TS, GVC_PB, E_sr}.
#' @export
bm_2023_trade_components <- function(io) {
  stopifnot(inherits(io, "bm_io"))

  bil <- bm_2023_tripartite_bilateral(io)

  agg <- stats::aggregate(
    cbind(GVC_sr, GVC_PF, GVC_TS, GVC_PB, E_sr) ~ exporter,
    data = bil,
    FUN  = sum
  )

  agg
}

#' BM_2023 trade-based GVC participation measures
#'
#' From exporter-level GVC components, computes:
#' \itemize{
#'   \item \code{share_GVC_trade = GVC_sr / E_sr}
#'   \item \code{share_PF_trade  = GVC_PF / GVC_sr}
#'   \item \code{share_TS_trade  = GVC_TS / GVC_sr}
#'   \item \code{share_PB_trade  = GVC_PB / GVC_sr}
#'   \item \code{forward_trade   = (GVC_PF - GVC_PB) / GVC_sr}
#' }
#'
#' @param io bm_io object
#'
#' @return data.frame with exporter-level trade-based participation indicators.
#' @export
bm_2023_trade_measures <- function(io) {
  stopifnot(inherits(io, "bm_io"))

  comp <- bm_2023_trade_components(io)

  comp$share_GVC_trade <- comp$GVC_sr / comp$E_sr
  comp$share_PF_trade  <- comp$GVC_PF / comp$GVC_sr
  comp$share_TS_trade  <- comp$GVC_TS / comp$GVC_sr
  comp$share_PB_trade  <- comp$GVC_PB / comp$GVC_sr
  comp$forward_trade   <- (comp$GVC_PF - comp$GVC_PB) / comp$GVC_sr

  comp
}
