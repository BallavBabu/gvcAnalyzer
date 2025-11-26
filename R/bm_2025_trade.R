# R/bm_2025_trade.R

# --- Internal Helpers ---

bm_2025_DAVAX_sr <- function(io, s, r, L_list) {
  s <- bm_country_id(io, s)
  r <- bm_country_id(io, r)
  N <- io$N

  if (s == r) stop("s and r must be different")

  A    <- io$A
  Y    <- io$Y
  vvec <- io$v

  idx_s <- bm_idx_country(io, s)
  idx_r <- bm_idx_country(io, r)

  v_s  <- vvec[idx_s]
  L_ss <- L_list[[s]]
  L_rr <- L_list[[r]]
  A_sr <- A[idx_s, idx_r, drop = FALSE]
  Y_sr <- Y[idx_s, r]
  Y_rr <- Y[idx_r, r]

  term   <- Y_sr + A_sr %*% (L_rr %*% Y_rr)
  result <- as.numeric(v_s %*% (L_ss %*% term))
  result
}

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

  term_final <- Y_sr + A_sr %*% (L_rr %*% Y_rr)
  q_s        <- L_ss %*% term_final

  PB_val <- 0
  for (t in seq_len(G)) {
    if (t == s) next
    idx_t <- bm_idx_country(io, t)
    A_ts  <- A[idx_t, idx_s, drop = FALSE]
    PB_val <- PB_val + as.numeric(u_N %*% (A_ts %*% q_s))
  }
  PB_val
}

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

  sum_e_r_to_others <- rep(0, N)
  for (j in seq_len(G)) {
    if (j == r) next
    sum_e_r_to_others <- sum_e_r_to_others + bm_get_e_sr(io, r, j)
  }

  term_two <- A_sr %*% (L_rr %*% sum_e_r_to_others)
  q_s      <- L_ss %*% term_two

  TS_val <- 0
  for (t in seq_len(G)) {
    if (t == s) next
    idx_t <- bm_idx_country(io, t)
    A_ts  <- A[idx_t, idx_s, drop = FALSE]
    TS_val <- TS_val + as.numeric(u_N %*% (A_ts %*% q_s))
  }
  TS_val
}

# --- Exported Functions ---

#' BM_2025 tripartite GVC trade decomposition for one pair (s,r)
#' @param io A \code{bm_io} object.
#' @param s Exporter country (name or index).
#' @param r Importer country (name or index).
#' @return Data frame for the pair (s,r).
#' @export
bm_2025_tripartite_trade <- function(io, s, r) {
  stopifnot(inherits(io, "bm_io"))
  s <- bm_country_id(io, s)
  r <- bm_country_id(io, r)
  if (s == r) stop("s and r must be different")

  L_list <- bm_L_list(io)

  e_sr_vec <- bm_get_e_sr(io, s, r)
  E_sr     <- sum(e_sr_vec)

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
#' @param io A \code{bm_io} object.
#' @return Data frame for all pairs.
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

#' BM_2025 exporter-level GVC trade totals
#' @param io A \code{bm_io} object.
#' @return Data frame of exporter totals.
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
#' @param io A \code{bm_io} object.
#' @return Data frame of trade-based indicators.
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
