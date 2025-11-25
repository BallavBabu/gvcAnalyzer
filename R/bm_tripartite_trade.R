# bm_tripartite_trade.R
# Tripartite decomposition of GVC-related trade (BM_2025, Section 3.1)

# Internal: resolve country index from name or numeric
bm_country_id <- function(io, country) {
  if (is.numeric(country)) {
    stopifnot(country %in% seq_len(io$G))
    return(as.integer(country))
  }
  if (is.character(country)) {
    idx <- match(country, io$countries)
    if (is.na(idx)) {
      stop("Country '", country, "' not found in io$countries.")
    }
    return(as.integer(idx))
  }
  stop("country must be numeric index or character name.")
}

# Internal: row-of-ones of length N
bm_u_N <- function(io) {
  matrix(1, nrow = 1, ncol = io$N)
}

# Internal: shorthand for country block
bm_block_io <- function(io, i, j) {
  bm_block(io$A, i, j, io$N)
}

#' DAVAX_sr: Domestic value added absorbed in final demand in r
#'
#' Implements the DAVAX component for exports from s to r in BM_2025
#' tripartite GVC trade decomposition.
#'
#' It corresponds (in your notation) to:
#' \eqn{DAVAX_{sr} = v_s L_{ss} ( y_{sr} + A_{sr} L_{rr} y_{rr} )}.
#'
#' @param io A \code{bm_io} object created by \code{bm_build_io()}.
#' @param s Exporting country (index or name).
#' @param r Importing country (index or name).
#'
#' @return Numeric scalar, DAVAX_sr.
#' @export
bm_davax_sr <- function(io, s, r) {
  s_id <- bm_country_id(io, s)
  r_id <- bm_country_id(io, r)
  if (s_id == r_id) stop("s and r must be different countries.")

  N  <- io$N
  v  <- io$v
  Y  <- io$Y
  Ls <- io$L_list[[s_id]]
  Lr <- io$L_list[[r_id]]

  idx_s <- bm_idx_country(s_id, N)
  idx_r <- bm_idx_country(r_id, N)

  v_s  <- v[idx_s]                             # length N
  A_sr <- bm_block_io(io, s_id, r_id)         # N x N
  y_sr <- Y[idx_s, r_id, drop = FALSE]        # N x 1
  y_rr <- Y[idx_r, r_id, drop = FALSE]        # N x 1

  term <- y_sr + A_sr %*% (Lr %*% y_rr)       # N x 1
  as.numeric(v_s %*% (Ls %*% term))
}

#' PB_sr: Pure-backward GVC trade (import content used in r's final demand)
#'
#' Pure-backward GVC trade is the import content of production in s
#' that serves r's final demand of s and r goods.
#'
#' In your notation:
#' - q_s(sr) = L_ss ( y_{sr} + A_{sr} L_{rr} y_{rr} )
#' - PB_sr = Σ_{t≠s} ( u_N A_{ts} q_s(sr) ).
#'
#' @param io A \code{bm_io} object.
#' @param s Exporting country (index or name).
#' @param r Importing country (index or name).
#'
#' @return Numeric scalar, PB_sr.
#' @export
bm_pb_sr <- function(io, s, r) {
  s_id <- bm_country_id(io, s)
  r_id <- bm_country_id(io, r)
  if (s_id == r_id) stop("s and r must be different countries.")

  G  <- io$G
  N  <- io$N
  Y  <- io$Y
  Ls <- io$L_list[[s_id]]
  Lr <- io$L_list[[r_id]]
  uN <- bm_u_N(io)

  idx_s <- bm_idx_country(s_id, N)
  idx_r <- bm_idx_country(r_id, N)

  A_sr <- bm_block_io(io, s_id, r_id)          # N x N
  y_sr <- Y[idx_s, r_id, drop = FALSE]         # N x 1
  y_rr <- Y[idx_r, r_id, drop = FALSE]         # N x 1

  term_final <- y_sr + A_sr %*% (Lr %*% y_rr)  # N x 1
  q_s <- Ls %*% term_final                     # N x 1

  PB_val <- 0
  for (t_id in seq_len(G)) {
    if (t_id == s_id) next
    A_ts <- bm_block(io$A, t_id, s_id, N)      # N x N
    PB_val <- PB_val + as.numeric(uN %*% (A_ts %*% q_s))
  }
  PB_val
}

#' TS_sr: Two-sided GVC trade for exports from s to r
#'
#' Two-sided GVC trade is the import content in s used to serve r's exports
#' to third countries (forward-plus-backward chains).
#'
#' In your notation:
#' - sum_e_r_to_others = Σ_{j≠r} e_{rj}
#' - demand in s: A_{sr} L_{rr} sum_e_r_to_others
#' - q_s(ts) = L_ss A_{sr} L_{rr} sum_e_r_to_others
#' - TS_sr = Σ_{t≠s} ( u_N A_{ts} q_s(ts) ).
#'
#' @param io A \code{bm_io} object.
#' @param s Exporting country (index or name).
#' @param r Importing country (index or name).
#'
#' @return Numeric scalar, TS_sr.
#' @export
bm_ts_sr <- function(io, s, r) {
  s_id <- bm_country_id(io, s)
  r_id <- bm_country_id(io, r)
  if (s_id == r_id) stop("s and r must be different countries.")

  G  <- io$G
  N  <- io$N
  Ls <- io$L_list[[s_id]]
  Lr <- io$L_list[[r_id]]
  uN <- bm_u_N(io)

  A_sr <- bm_block_io(io, s_id, r_id)          # N x N

  # Σ_{j≠r} e_{rj}
  sum_e_r_to_others <- rep(0, N)
  for (j_id in seq_len(G)) {
    if (j_id == r_id) next
    sum_e_r_to_others <- sum_e_r_to_others + bm_get_e_sr(io, r_id, j_id)
  }

  term_two <- A_sr %*% (Lr %*% sum_e_r_to_others)  # N x 1
  q_s <- Ls %*% term_two                           # N x 1

  TS_val <- 0
  for (t_id in seq_len(G)) {
    if (t_id == s_id) next
    A_ts <- bm_block(io$A, t_id, s_id, N)          # N x N
    TS_val <- TS_val + as.numeric(uN %*% (A_ts %*% q_s))
  }
  TS_val
}

#' Tripartite GVC trade decomposition for a bilateral pair (s,r)
#'
#' For exports from s to r, this function returns:
#' \itemize{
#'   \item E_sr: gross exports (sum of e_sr)
#'   \item DAVAX_sr: domestic value added that crosses only one border and is absorbed in r's final demand
#'   \item GVC_sr: GVC-related part of E_sr, GVC_sr = E_sr - DAVAX_sr
#'   \item GVC_PF: pure-forward GVC trade (residual)
#'   \item GVC_TS: two-sided GVC trade
#'   \item GVC_PB: pure-backward GVC trade
#' }
#'
#' In your notation:
#' \deqn{E_{sr} = DAVAX_{sr} + GVC_{sr}}
#' \deqn{GVC_{sr} = GVC^{PF}_{sr} + GVC^{TS}_{sr} + GVC^{PB}_{sr}}.
#'
#' @param io A \code{bm_io} object.
#' @param s Exporting country (index or name).
#' @param r Importing country (index or name).
#'
#' @return A one-row \code{data.frame} with all components.
#' @export
bm_tripartite_trade_sr <- function(io, s, r) {
  s_id <- bm_country_id(io, s)
  r_id <- bm_country_id(io, r)
  if (s_id == r_id) stop("s and r must be different countries.")

  e_sr_vec <- bm_get_e_sr(io, s_id, r_id)
  E_sr     <- sum(e_sr_vec)

  DAVAX_sr <- bm_davax_sr(io, s_id, r_id)
  PB_sr    <- bm_pb_sr(io, s_id, r_id)
  TS_sr    <- bm_ts_sr(io, s_id, r_id)

  GVC_sr   <- E_sr - DAVAX_sr
  GVC_PB   <- PB_sr
  GVC_TS   <- TS_sr
  GVC_PF   <- GVC_sr - GVC_PB - GVC_TS

  data.frame(
    exporter = io$countries[s_id],
    importer = io$countries[r_id],
    E_sr     = E_sr,
    DAVAX_sr = DAVAX_sr,
    GVC_sr   = GVC_sr,
    GVC_PF   = GVC_PF,
    GVC_TS   = GVC_TS,
    GVC_PB   = GVC_PB
  )
}

#' Tripartite GVC trade for all bilateral pairs
#'
#' Computes the BM_2025 tripartite GVC trade decomposition for all
#' ordered pairs (s,r) with s ≠ r.
#'
#' @param io A \code{bm_io} object.
#'
#' @return A \code{data.frame} with one row per pair (s,r).
#' @export
bm_tripartite_trade_all <- function(io) {
  G <- io$G
  out_list <- list()
  k <- 1L
  for (s_id in seq_len(G)) {
    for (r_id in seq_len(G)) {
      if (s_id == r_id) next
      out_list[[k]] <- bm_tripartite_trade_sr(io, s_id, r_id)
      k <- k + 1L
    }
  }
  do.call(rbind, out_list)
}
