# R/bm_2025_trade.R
#
# BM_2025 tripartite decomposition of GVC-related trade (PF/TS/PB)
# and trade-based participation indicators.
#
# This implements the 3.1 + 4.1 logic:
#   E_sr = DAVAX_sr + GVC_sr
#   GVC_sr = GVC_PF_sr + GVC_TS_sr + GVC_PB_sr

# Internal short-hands for indices
bm_idx_country <- function(io, g) {
  io$.idx_country(g)
}

# Internal: build domestic Leontief inverses L_gg = (I - A_gg)^{-1}
bm_L_list <- function(io) {
  G <- io$G
  N <- io$N
  A <- io$A
  L_list <- vector("list", G)
  I_N <- Matrix::Diagonal(n = N)
  for (g in seq_len(G)) {
    idx_g <- bm_idx_country(io, g)
    A_gg  <- A[idx_g, idx_g, drop = FALSE]
    L_list[[g]] <- solve(I_N - A_gg)
  }
  L_list
}

# -------------------------------------------------------------------
#  BM_2025 building blocks for one bilateral pair (s,r)
# -------------------------------------------------------------------

# Internal: DAVAX_sr (domestic VA of s that crosses one border and is
# absorbed in r's final demand of s- and r-produced goods).
bm_2025_DAVAX_sr <- function(io, s, r, L_list) {
  s <- bm_country_id(io, s)
  r <- bm_country_id(io, r)
  G <- io$G
  N <- io$N

  if (s == r) stop("s and r must be different")

  A    <- io$A
  Y    <- io$Y
  vvec <- io$v

  idx_s <- bm_idx_country(io, s)
  idx_r <- bm_idx_country(io, r)

  v_s  <- vvec[idx_s]                 # length N
  L_ss <- L_list[[s]]                 # N x N
  L_rr <- L_list[[r]]                 # N x N
  A_sr <- A[idx_s, idx_r, drop = FALSE]
  Y_sr <- Y[idx_s, r]                 # N x 1
  Y_rr <- Y[idx_r, r]                 # N x 1

  term   <- Y_sr + A_sr %*% (L_rr %*% Y_rr)  # N x 1
  result <- as.numeric(v_s %*% (L_ss %*% term))
  result
}

# Internal: pure-backward GVC trade PB_sr
#   = gross import content of exports from s->r that stop in r's final demand.
bm_2025_PB_sr <- function(io, s, r, L_list) {
  s <- bm_country_id(io, s)
  r <- bm_country_id(io, r)
  G <- io$G
  N <- io$N

  if (s == r) stop("s and r must be different")

  A <- io$A
  Y <- io$Y
  u_N <- matrix(1, nrow = 1, ncol = N)

  idx_s <- bm_idx_country(io, s)
  idx_r <- bm_idx_country(io, r)

  L_ss <- L_list[[s]]
  L_rr <- L_list[[r]]
  A_sr <- A[idx_s, idx_r, drop = FALSE]
  Y_sr <- Y[idx_s, r]
  Y_rr <- Y[idx_r, r]

  term_final <- Y_sr + A_sr %*% (L_rr %*% Y_rr)   # N x 1
  q_s        <- L_ss %*% term_final               # N x 1

  PB_val <- 0
  for (t in seq_len(G)) {
    if (t == s) next
    idx_t <- bm_idx_country(io, t)
    A_ts  <- A[idx_t, idx_s, drop = FALSE]
    PB_val <- PB_val + as.numeric(u_N %*% (A_ts %*% q_s))
  }
  PB_val
}

# Internal: two-sided GVC trade TS_sr
#   = imported inputs in s used to produce exports from s->r,
#     which r uses to export to other partners.
bm_2025_TS_sr <- function(io, s, r, L_list) {
  s <- bm_country_id(io, s)
  r <- bm_country_id(io, r)
  G <- io$G
  N <- io$N

  if (s == r) stop("s and r must be different")

  A <- io$A
  Y <- io$Y
  u_N <- matrix(1, nrow = 1, ncol = N)

  idx_s <- bm_idx_country(io, s)
  idx_r <- bm_idx_country(io, r)

  L_ss <- L_list[[s]]
  L_rr <- L_list[[r]]
  A_sr <- A[idx_s, idx_r, drop = FALSE]

  # Exports of r to all j ≠ r (intermediate + final)
  sum_e_r_to_others <- rep(0, N)
  for (j in seq_len(G)) {
    if (j == r) next
    sum_e_r_to_others <- sum_e_r_to_others + bm_e_sr(io, r, j)
  }

  term_two <- A_sr %*% (L_rr %*% sum_e_r_to_others)  # N x 1
  q_s      <- L_ss %*% term_two                      # N x 1

  TS_val <- 0
  for (t in seq_len(G)) {
    if (t == s) next
    idx_t <- bm_idx_country(io, t)
    A_ts  <- A[idx_t, idx_s, drop = FALSE]
    TS_val <- TS_val + as.numeric(u_N %*% (A_ts %*% q_s))
  }
  TS_val
}

# -------------------------------------------------------------------
#  Public: BM_2025 bilateral tripartite decomposition (trade)
# -------------------------------------------------------------------

#' BM_2025 tripartite GVC trade decomposition for one pair (s,r)
#'
#' Decompose gross exports \eqn{E_{sr}} (from exporter s to importer r)
#' into:
#'
#'   E_{sr} = DAVAX_{sr} + GVC_{sr},
#'   GVC_{sr} = GVC_PF_{sr} + GVC_TS_{sr} + GVC_PB_{sr}.
#'
#' This follows the BM_2025 tripartite concept:
#'   \itemize{
#'     \item DAVAX_sr: domestic VA of s crossing one border and absorbed
#'           in r's final demand (of s and r goods).
#'     \item GVC_PB_sr: pure backward GVC trade (import content of
#'           exports from s->r that stop in r's final demand).
#'     \item GVC_TS_sr: two-sided GVC trade (imported inputs in s used to
#'           produce exports from s->r that r uses to export further).
#'     \item GVC_PF_sr: pure forward GVC trade, as residual:
#'           GVC_PF_sr = GVC_sr - GVC_TS_sr - GVC_PB_sr.
#'   }
#'
#' @param io bm_io object
#' @param s  exporter (country index or code, e.g. 1 or "China")
#' @param r  importer (country index or code, e.g. 2 or "India")
#'
#' @return A data frame with one row for the pair (s,r):
#'   exporter, importer, E_sr, DAVAX_sr, GVC_sr, GVC_PF, GVC_TS, GVC_PB.
#' @export
bm_2025_tripartite_trade <- function(io, s, r) {
  stopifnot(inherits(io, "bm_io"))
  s <- bm_country_id(io, s)
  r <- bm_country_id(io, r)
  if (s == r) stop("s and r must be different")

  G <- io$G

  # Precompute domestic Leontief once (for now, inside the function)
  L_list <- bm_L_list(io)

  # Gross exports
  e_sr_vec <- bm_e_sr(io, s, r)
  E_sr     <- sum(e_sr_vec)

  # DAVAX, PB, TS
  DAVAX_sr <- bm_2025_DAVAX_sr(io, s, r, L_list)
  PB_sr    <- bm_2025_PB_sr(io, s, r, L_list)
  TS_sr    <- bm_2025_TS_sr(io, s, r, L_list)

  GVC_sr <- E_sr - DAVAX_sr
  GVC_PF <- GVC_sr - PB_sr - TS_sr

  data.frame(
    exporter = io$countries[s],
    importer = io$countries[r],
    E_sr     = E_sr,
    DAVAX_sr = DAVAX_sr,
    GVC_sr   = GVC_sr,
    GVC_PF   = GVC_PF,
    GVC_TS   = TS_sr,
    GVC_PB   = PB_sr
  )
}

#' BM_2025 tripartite GVC trade decomposition for all pairs
#'
#' Applies \code{bm_2025_tripartite_trade()} to all ordered pairs
#' (s,r), s ≠ r.
#'
#' @param io bm_io object
#'
#' @return data.frame with one row per exporter-importer pair:
#'   exporter, importer, E_sr, DAVAX_sr, GVC_sr, GVC_PF, GVC_TS, GVC_PB.
#' @export
bm_2025_tripartite_trade_all <- function(io) {
  stopifnot(inherits(io, "bm_io"))
  G <- io$G

  records <- list()
  k <- 1L

  for (s in seq_len(G)) {
    for (r in seq_len(G)) {
      if (s == r) next
      records[[k]] <- bm_2025_tripartite_trade(io, s, r)
      k <- k + 1L
    }
  }

  out <- do.call(rbind, records)
  rownames(out) <- NULL
  out
}

# -------------------------------------------------------------------
#  Exporter-level aggregation and trade participation indicators
# -------------------------------------------------------------------

#' BM_2025 exporter-level GVC trade totals
#'
#' Aggregates BM_2025 tripartite trade decomposition over all destinations
#' for each exporter s:
#'
#'   \eqn{E_s = sum_r E_{sr}},
#'   \eqn{GVC_s = sum_r GVC_{sr}},
#'   and similarly for the PF/TS/PB components.
#'
#' @param io bm_io object
#'
#' @return data.frame with one row per exporter:
#'   exporter, E_s, GVC_s, GVC_PF_s, GVC_TS_s, GVC_PB_s.
#' @export
bm_2025_trade_exporter <- function(io) {
  bilateral <- bm_2025_tripartite_trade_all(io)
  agg <- stats::aggregate(
    cbind(E_sr, GVC_sr, GVC_PF, GVC_TS, GVC_PB) ~ exporter,
    data = bilateral,
    FUN  = sum
  )
  names(agg) <- c("exporter", "E_s", "GVC_s", "GVC_PF_s", "GVC_TS_s", "GVC_PB_s")
  agg
}

#' BM_2025 trade-based GVC participation indicators
#'
#' Computes GVC participation and composition indicators for each exporter:
#'
#'   \itemize{
#'     \item share_GVC_trade = GVC_s / E_s
#'     \item share_PF_trade  = GVC_PF_s / GVC_s
#'     \item share_TS_trade  = GVC_TS_s / GVC_s
#'     \item share_PB_trade  = GVC_PB_s / GVC_s
#'     \item forward_trade   = (GVC_PF_s - GVC_PB_s) / GVC_s
#'   }
#'
#' @param io bm_io object
#'
#' @return data.frame with exporter-level totals and indicator columns.
#' @export
bm_2025_trade_measures <- function(io) {
  agg <- bm_2025_trade_exporter(io)

  agg$share_GVC_trade <- agg$GVC_s / agg$E_s
  agg$share_PF_trade  <- agg$GVC_PF_s / agg$GVC_s
  agg$share_TS_trade  <- agg$GVC_TS_s / agg$GVC_s
  agg$share_PB_trade  <- agg$GVC_PB_s / agg$GVC_s
  agg$forward_trade   <- (agg$GVC_PF_s - agg$GVC_PB_s) / agg$GVC_s

  agg
}
