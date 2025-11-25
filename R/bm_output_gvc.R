# bm_output_gvc.R
# Output-based GVC concepts and measures (BM_2025 Section 3.2)

# We reuse:
# - bm_country_id(io, country)    ## defined in bm_tripartite_trade.R
# - bm_idx_country(g, N)         ## defined in bm_io.R
# - bm_block(M, i, j)         ## defined in bm_io.R
# - bm_get_e_star(io, s)         ## from bm_io.R

# --- Small helpers -------------------------------------------------------

bm_v_country <- function(io, g) {
  g_id <- bm_country_id(io, g)
  idx  <- bm_idx_country(g_id, io$N)
  io$v[idx]
}

bm_X_country <- function(io, g) {
  g_id <- bm_country_id(io, g)
  idx  <- bm_idx_country(g_id, io$N)
  io$X[idx]
}

bm_Y_country_prod_to_dest <- function(io, g, dest) {
  g_id   <- bm_country_id(io, g)
  dest_id <- bm_country_id(io, dest)
  idx_g  <- bm_idx_country(g_id, io$N)
  io$Y[idx_g, dest_id, drop = FALSE]   # N x 1
}

bm_Y_country_total <- function(io, g) {
  g_id  <- bm_country_id(io, g)
  idx_g <- bm_idx_country(g_id, io$N)
  rowSums(io$Y[idx_g, , drop = FALSE]) # N x 1
}

# Precompute export-related output Xexp_r = B_rr e_r*
bm_build_Xexp_list <- function(io) {
  G <- io$G; N <- io$N
  Xexp_list <- vector("list", G)
  for (r_id in seq_len(G)) {
    B_rr      <- bm_block(io$B, r_id, r_id)
    e_r_star  <- bm_get_e_star(io, r_id)
    Xexp_list[[r_id]] <- B_rr %*% e_r_star
  }
  Xexp_list
}

# --- 3.2.1 Pure-forward GVC-related output of s ---------------------------

#' Pure-forward GVC-related output (country level)
#'
#' For a given country s, this function computes the pure-forward
#' GVC-related output, i.e. domestic value added of s embodied in r's
#' export-related output to third countries (forward chains).
#'
#' @param io A \code{bm_io} object.
#' @param s  Country (index or name).
#' @param Xexp_list Optional precomputed list of export-related output
#'   (as from \code{bm_build_Xexp_list(io)}). If \code{NULL}, it is
#'   built internally.
#'
#' @return Numeric scalar: pure-forward GVC-related output of s.
#' @export
bm_gvc_pf_output <- function(io, s, Xexp_list = NULL) {
  s_id <- bm_country_id(io, s)
  G    <- io$G
  N    <- io$N

  if (is.null(Xexp_list)) {
    Xexp_list <- bm_build_Xexp_list(io)
  }

  v_s  <- matrix(bm_v_country(io, s_id), nrow = 1)
  L_ss <- io$L_list[[s_id]]
  A_ss <- bm_block(io$A, s_id, s_id)

  val <- 0
  for (r_id in seq_len(G)) {
    if (r_id == s_id) next
    A_sr   <- bm_block(io$A, s_id, r_id)
    Xexp_r <- Xexp_list[[r_id]]               # N x 1

    inner   <- A_sr %*% Xexp_r                # direct sales from s to r's export output
    inner_u <- A_ss %*% (L_ss %*% inner)      # upstream domestic propagation

    q_sr <- inner + inner_u                   # N x 1 gross output in s
    val  <- val + as.numeric(v_s %*% (L_ss %*% q_sr))
  }
  val
}

# --- 3.2.2 Pure-backward GVC-related output of s --------------------------

#' Pure-backward GVC-related output (country level)
#'
#' For a given country s, this function computes the pure-backward
#' GVC-related output, interpreted as foreign value added embodied
#' in s's output that is ultimately absorbed abroad in a one-border
#' GVC pattern.
#'
#' This follows the BM_2025 style expression using:
#' - Y_s_total = Σ_z Y_sz
#' - Y_ss      = Y_ss
#' - A_jk, B_ks, L_jj, L_ss blocks.
#'
#' @param io A \code{bm_io} object.
#' @param s  Country (index or name).
#'
#' @return Numeric scalar: pure-backward GVC-related output of s.
#' @export
bm_gvc_pb_output <- function(io, s) {
  s_id <- bm_country_id(io, s)
  G    <- io$G
  N    <- io$N

  Y_s_total <- bm_Y_country_total(io, s_id)            # N x 1
  Y_ss      <- bm_Y_country_prod_to_dest(io, s_id, s_id)

  val_first  <- 0
  val_second <- 0

  for (j_id in seq_len(G)) {
    v_j  <- matrix(bm_v_country(io, j_id), nrow = 1)
    L_jj <- io$L_list[[j_id]]

    # Σ_{k≠j} A_jk B_ks Y_s_total
    term_j <- rep(0, N)
    for (k_id in seq_len(G)) {
      if (k_id == j_id) next
      A_jk <- bm_block(io$A, j_id, k_id)
      B_ks <- bm_block(io$B, k_id, s_id)
      term_j <- term_j + A_jk %*% (B_ks %*% Y_s_total)
    }
    val_first <- val_first + as.numeric(v_j %*% (L_jj %*% term_j))

    # subtract domestic chain part for j ≠ s
    if (j_id != s_id) {
      A_js <- bm_block(io$A, j_id, s_id)
      L_ss <- io$L_list[[s_id]]
      val_second <- val_second + as.numeric(
        v_j %*% (L_jj %*% (A_js %*% (L_ss %*% Y_ss)))
      )
    }
  }
  val_first - val_second
}

# --- 3.2.3 Two-sided GVC-related output of s (imported & domestic) -------

#' Two-sided GVC-related output from imported inputs
#'
#' Imported two-sided GVC-related output of s is defined as all imported
#' inputs used in s's output, net of those absorbed in purely domestic
#' chains, net of the pure-backward part.
#'
#' @param io A \code{bm_io} object.
#' @param s  Country (index or name).
#' @param gvc_pb_s Pre-computed pure-backward GVC output of s
#'   (from \code{bm_gvc_pb_output}).
#'
#' @return Numeric scalar: imported two-sided GVC-related output of s.
#' @export
bm_gvc_ts_import_output <- function(io, s, gvc_pb_s) {
  s_id <- bm_country_id(io, s)
  G    <- io$G
  N    <- io$N

  X_s  <- bm_X_country(io, s_id)
  Y_ss <- bm_Y_country_prod_to_dest(io, s_id, s_id)

  val_first  <- 0
  val_second <- 0

  for (j_id in seq_len(G)) {
    v_j  <- matrix(bm_v_country(io, j_id), nrow = 1)
    L_jj <- io$L_list[[j_id]]

    # Σ_{k≠j} A_jk B_ks X_s
    term_j1 <- rep(0, N)
    for (k_id in seq_len(G)) {
      if (k_id == j_id) next
      A_jk <- bm_block(io$A, j_id, k_id)
      B_ks <- bm_block(io$B, k_id, s_id)
      term_j1 <- term_j1 + A_jk %*% (B_ks %*% X_s)
    }
    val_first <- val_first + as.numeric(v_j %*% (L_jj %*% term_j1))

    # subtract purely domestic-chain part for j ≠ s
    if (j_id != s_id) {
      A_js <- bm_block(io$A, j_id, s_id)
      L_ss <- io$L_list[[s_id]]
      val_second <- val_second + as.numeric(
        v_j %*% (L_jj %*% (A_js %*% (L_ss %*% (L_ss %*% Y_ss))))
      )
    }
  }
  # Imported two-sided = all imported inputs in X_s
  # minus those in pure domestic chains and pure-backward
  val_first - val_second - gvc_pb_s
}

#' Two-sided GVC-related output from domestic inputs
#'
#' Domestic two-sided GVC-related output of s is the domestic part of
#' two-sided GVC chains: domestic value added of s embodied in r's GVC
#' production that combines domestic and foreign inputs.
#'
#' @param io A \code{bm_io} object.
#' @param s  Country (index or name).
#' @param Xexp_list Optional precomputed export-related outputs list.
#'
#' @return Numeric scalar: domestic two-sided GVC-related output of s.
#' @export
bm_gvc_ts_domestic_output <- function(io, s, Xexp_list = NULL) {
  s_id <- bm_country_id(io, s)
  G    <- io$G
  N    <- io$N

  if (is.null(Xexp_list)) {
    Xexp_list <- bm_build_Xexp_list(io)
  }

  v_s  <- matrix(bm_v_country(io, s_id), nrow = 1)
  L_ss <- io$L_list[[s_id]]
  A_ss <- bm_block(io$A, s_id, s_id)

  val <- 0
  for (r_id in seq_len(G)) {
    if (r_id == s_id) next
    A_sr   <- bm_block(io$A, s_id, r_id)
    Xexp_r <- Xexp_list[[r_id]]

    inner  <- A_sr %*% Xexp_r
    inner2 <- A_ss %*% (L_ss %*% inner)
    q_sr   <- inner + inner2
    val    <- val + as.numeric(v_s %*% (L_ss %*% (A_ss %*% q_sr)))
  }
  val
}

# --- Country-level wrapper: all output-based components -------------------

#' Output-based GVC decomposition for one country
#'
#' For a given country s, this function returns a row with:
#' \itemize{
#'   \item GVC_PF_X  : pure-forward GVC-related output
#'   \item GVC_PB_X  : pure-backward GVC-related output
#'   \item GVC_TSImp : imported-input two-sided GVC-related output
#'   \item GVC_TSDom : domestic-input two-sided GVC-related output
#'   \item GVC_TS_X  : total two-sided (sum of imported and domestic)
#'   \item GVC_X     : total GVC-related output (PF + PB + TS)
#'   \item DomX      : proxy for purely domestic output chains
#'   \item TradX     : residual "traditional" (one-border) trade-related output
#'   \item X_total   : total output of s
#' }
#'
#' @param io A \code{bm_io} object.
#' @param s  Country (index or name).
#' @param Xexp_list Optional precomputed export-related outputs list.
#'
#' @return A one-row \code{data.frame} with these components.
#' @export
bm_gvc_output_country <- function(io, s, Xexp_list = NULL) {
  s_id <- bm_country_id(io, s)
  N    <- io$N

  if (is.null(Xexp_list)) {
    Xexp_list <- bm_build_Xexp_list(io)
  }

  PF_X_s <- bm_gvc_pf_output(io, s_id, Xexp_list)
  PB_X_s <- bm_gvc_pb_output(io, s_id)
  TS_imp <- bm_gvc_ts_import_output(io, s_id, PB_X_s)
  TS_dom <- bm_gvc_ts_domestic_output(io, s_id, Xexp_list)
  TS_X_s <- TS_imp + TS_dom
  GVC_X_s <- PF_X_s + PB_X_s + TS_X_s

  X_s_vec  <- bm_X_country(io, s_id)
  X_s_tot  <- sum(X_s_vec)

  # Simple, transparent proxy: domestic output as domestic final demand output
  Y_ss  <- bm_Y_country_prod_to_dest(io, s_id, s_id)
  v_s   <- matrix(bm_v_country(io, s_id), nrow = 1)
  L_ss  <- io$L_list[[s_id]]
  DomX  <- as.numeric(v_s %*% (L_ss %*% Y_ss))

  TradX <- X_s_tot - DomX - GVC_X_s

  data.frame(
    country   = io$countries[s_id],
    GVC_PF_X  = PF_X_s,
    GVC_PB_X  = PB_X_s,
    GVC_TSImp = TS_imp,
    GVC_TSDom = TS_dom,
    GVC_TS_X  = TS_X_s,
    GVC_X     = GVC_X_s,
    DomX      = DomX,
    TradX     = TradX,
    X_total   = X_s_tot
  )
}

#' Output-based GVC decomposition for all countries
#'
#' @param io A \code{bm_io} object.
#'
#' @return A \code{data.frame} with one row per country.
#' @export
bm_gvc_output_all <- function(io) {
  G          <- io$G
  Xexp_list  <- bm_build_Xexp_list(io)
  out_list   <- vector("list", G)
  for (s_id in seq_len(G)) {
    out_list[[s_id]] <- bm_gvc_output_country(io, s_id, Xexp_list)
  }
  do.call(rbind, out_list)
}
