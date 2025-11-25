# R/bm_2025_output.R
#
# BM_2025 output-based GVC measures (PF/TS/PB) and indicators.
# Self-contained: does not rely on bm_idx_country() or bm_L_list().

# ------------------------------------------------------------
# Local helpers: indices and country-level slices
# ------------------------------------------------------------

bm25_idx_country <- function(io, g) {
  # g can be numeric index or country name
  if (is.character(g)) {
    g <- match(g, io$countries)
  }
  if (is.na(g) || g < 1L || g > io$G) {
    stop("Invalid country identifier 'g' in bm25_idx_country().")
  }
  ((g - 1L) * io$N + 1L):(g * io$N)
}

bm25_v_country <- function(io, g) {
  idx <- bm25_idx_country(io, g)
  io$v[idx]
}

bm25_X_country <- function(io, g) {
  idx <- bm25_idx_country(io, g)
  io$X[idx]
}

bm25_Y_country_prod_to_dest <- function(io, g, dest) {
  idx <- bm25_idx_country(io, g)
  io$Y[idx, dest]
}

bm25_Y_country_total <- function(io, g) {
  idx <- bm25_idx_country(io, g)
  rowSums(io$Y[idx, , drop = FALSE])
}

# Domestic Leontief matrices L_gg
bm25_L_list <- function(io) {
  G <- io$G
  N <- io$N
  A <- io$A

  L_list <- vector("list", G)
  for (g in seq_len(G)) {
    idx_g <- bm25_idx_country(io, g)
    A_gg  <- as.matrix(A[idx_g, idx_g, drop = FALSE])
    L_list[[g]] <- solve(Matrix::Diagonal(N) - A_gg)
  }
  L_list
}

# -------------------------------------------------------------------
# Export-related output Xexp_r for each country r
#   e_r* = intermediate exports + final exports from r to all z ≠ r
#   Xexp_r = B_rr * e_r*
# -------------------------------------------------------------------
bm_2025_Xexp_list <- function(io) {
  G     <- io$G
  N     <- io$N
  GN    <- io$GN
  Z     <- io$Z
  Y     <- io$Y
  B_all <- io$B

  Xexp_list <- vector("list", G)

  for (r in seq_len(G)) {
    idx_r  <- bm25_idx_country(io, r)

    # e_r* (N-vector)
    rows_r  <- idx_r
    int_exp <- rowSums(
      Z[rows_r, setdiff(seq_len(GN), idx_r), drop = FALSE]
    )
    fin_exp <- rowSums(
      Y[rows_r, setdiff(seq_len(G), r), drop = FALSE]
    )
    e_rstar <- as.numeric(int_exp + fin_exp)   # length N

    # B_rr (N x N)
    B_rr <- as.matrix(B_all[idx_r, idx_r, drop = FALSE])

    # Xexp_r (N x 1)
    Xexp_list[[r]] <- as.numeric(B_rr %*% e_rstar)
  }

  Xexp_list
}

# ------------------------------------------------------------
# 3.2 Output-based PF / PB / TS components for one country s
# ------------------------------------------------------------

# 3.2.1 Pure-forward GVC-related output of s
bm_2025_GVC_PF_X <- function(io, s, L_list, Xexp_list) {
  if (is.character(s)) s <- match(s, io$countries)
  G <- io$G
  A <- io$A

  idx_s <- bm25_idx_country(io, s)

  v_s  <- matrix(bm25_v_country(io, s), nrow = 1)  # 1 x N
  L_ss <- L_list[[s]]                              # N x N
  A_ss <- as.matrix(A[idx_s, idx_s, drop = FALSE])

  val <- 0
  for (r in seq_len(G)) {
    if (r == s) next
    idx_r   <- bm25_idx_country(io, r)
    A_sr    <- as.matrix(A[idx_s, idx_r, drop = FALSE])     # s->r
    Xexp_r  <- matrix(Xexp_list[[r]], ncol = 1)             # N x 1

    inner   <- A_sr %*% Xexp_r                              # N x 1
    inner_u <- A_ss %*% (L_ss %*% inner)                    # N x 1

    q_sr <- inner + inner_u                                 # N x 1
    val  <- val + as.numeric(v_s %*% (L_ss %*% q_sr))
  }
  val
}

# 3.2.2 Pure-backward GVC-related output of s
bm_2025_GVC_PB_X <- function(io, s, L_list) {
  if (is.character(s)) s <- match(s, io$countries)
  G <- io$G
  A <- io$A
  B <- io$B

  Y_s_total <- matrix(bm25_Y_country_total(io, s), ncol = 1)       # N x 1
  Y_ss      <- matrix(bm25_Y_country_prod_to_dest(io, s, s), ncol = 1)

  val_first  <- 0
  val_second <- 0

  for (j in seq_len(G)) {
    idx_j  <- bm25_idx_country(io, j)
    v_j    <- matrix(bm25_v_country(io, j), nrow = 1)
    L_jj   <- L_list[[j]]

    # Σ_{k≠j} A_jk B_ks Y_s_total
    term_j <- matrix(0, nrow = length(idx_j), ncol = 1)
    for (k in seq_len(G)) {
      if (k == j) next
      idx_k <- bm25_idx_country(io, k)
      A_jk  <- as.matrix(A[idx_j, idx_k, drop = FALSE])
      B_ks  <- as.matrix(B[idx_k, bm25_idx_country(io, s), drop = FALSE])
      term_j <- term_j + A_jk %*% (B_ks %*% Y_s_total)
    }
    val_first <- val_first + as.numeric(v_j %*% (L_jj %*% term_j))

    # second term: j ≠ s
    if (j != s) {
      A_js <- as.matrix(A[idx_j, bm25_idx_country(io, s), drop = FALSE])
      L_ss <- L_list[[s]]
      val_second <- val_second +
        as.numeric(v_j %*% (L_jj %*% (A_js %*% (L_ss %*% Y_ss))))
    }
  }

  val_first - val_second
}

# 3.2.3 Two-sided GVC-related output of s, imported-input component
bm_2025_GVC_TS_Imp_X <- function(io, s, L_list, GVC_PB_X_s) {
  if (is.character(s)) s <- match(s, io$countries)
  G <- io$G
  A <- io$A
  B <- io$B

  X_s  <- matrix(bm25_X_country(io, s), ncol = 1)
  Y_ss <- matrix(bm25_Y_country_prod_to_dest(io, s, s), ncol = 1)

  val_first  <- 0
  val_second <- 0

  for (j in seq_len(G)) {
    idx_j  <- bm25_idx_country(io, j)
    v_j    <- matrix(bm25_v_country(io, j), nrow = 1)
    L_jj   <- L_list[[j]]

    # Σ_{k≠j} A_jk B_ks X_s
    term_j1 <- matrix(0, nrow = length(idx_j), ncol = 1)
    for (k in seq_len(G)) {
      if (k == j) next
      idx_k <- bm25_idx_country(io, k)
      A_jk  <- as.matrix(A[idx_j, idx_k, drop = FALSE])
      B_ks  <- as.matrix(B[idx_k, bm25_idx_country(io, s), drop = FALSE])
      term_j1 <- term_j1 + A_jk %*% (B_ks %*% X_s)
    }
    val_first <- val_first + as.numeric(v_j %*% (L_jj %*% term_j1))

    # second term: j ≠ s
    if (j != s) {
      A_js <- as.matrix(A[idx_j, bm25_idx_country(io, s), drop = FALSE])
      L_ss <- L_list[[s]]
      val_second <- val_second +
        as.numeric(
          v_j %*% (L_jj %*% (A_js %*% (L_ss %*% (L_ss %*% Y_ss))))
        )
    }
  }

  # Imported-input two-sided = all imported inputs in X_s
  # minus those absorbed in "pure domestic" chains, minus pure-backward
  val_first - val_second - GVC_PB_X_s
}

# 3.2.4 Two-sided GVC-related output of s, domestic-input component
bm_2025_GVC_TS_Dom_X <- function(io, s, L_list, Xexp_list) {
  if (is.character(s)) s <- match(s, io$countries)
  G <- io$G
  A <- io$A

  idx_s <- bm25_idx_country(io, s)

  v_s  <- matrix(bm25_v_country(io, s), nrow = 1)
  L_ss <- L_list[[s]]
  A_ss <- as.matrix(A[idx_s, idx_s, drop = FALSE])

  val <- 0
  for (r in seq_len(G)) {
    if (r == s) next
    idx_r   <- bm25_idx_country(io, r)
    A_sr    <- as.matrix(A[idx_s, idx_r, drop = FALSE])
    Xexp_r  <- matrix(Xexp_list[[r]], ncol = 1)
    inner   <- A_sr %*% Xexp_r
    inner2  <- A_ss %*% (L_ss %*% inner)
    q_sr    <- inner + inner2
    val     <- val + as.numeric(v_s %*% (L_ss %*% (A_ss %*% q_sr)))
  }
  val
}

# ------------------------------------------------------------
# Aggregate for each exporter s
# ------------------------------------------------------------

#' BM_2025 output-based GVC components by exporter
#'
#' Computes output-based PF/TS/PB GVC components for each country s.
#'
#' @param io bm_io object
#'
#' @return data.frame with one row per country.
#' @export
bm_2025_output_components <- function(io) {
  stopifnot(inherits(io, "bm_io"))
  G <- io$G

  # Domestic Leontief L_gg and export-related output in each r
  L_list    <- bm25_L_list(io)
  Xexp_list <- bm_2025_Xexp_list(io)

  output_list <- vector("list", G)

  for (s in seq_len(G)) {
    # core GVC components
    PF_X_s <- bm_2025_GVC_PF_X(io, s, L_list, Xexp_list)
    PB_X_s <- bm_2025_GVC_PB_X(io, s, L_list)
    TS_imp <- bm_2025_GVC_TS_Imp_X(io, s, L_list, PB_X_s)
    TS_dom <- bm_2025_GVC_TS_Dom_X(io, s, L_list, Xexp_list)
    TS_X_s <- TS_imp + TS_dom
    GVC_X_s <- PF_X_s + PB_X_s + TS_X_s

    # total output of s
    X_s_tot <- sum(bm25_X_country(io, s))

    # simple domestic-output proxy (BM-style toy)
    idx_s <- bm25_idx_country(io, s)
    Y_ss  <- matrix(bm25_Y_country_prod_to_dest(io, s, s), ncol = 1)
    v_s   <- matrix(bm25_v_country(io, s), nrow = 1)
    L_ss  <- L_list[[s]]
    A_ss  <- as.matrix(io$A[idx_s, idx_s, drop = FALSE])

    DomX_s <- as.numeric(
      v_s %*% (L_ss %*% Y_ss) +
        v_s %*% (L_ss %*% (A_ss %*% (L_ss %*% Y_ss)))
    )

    TradX_s <- X_s_tot - DomX_s - GVC_X_s

    output_list[[s]] <- data.frame(
      country   = io$countries[s],
      GVC_PF_X  = PF_X_s,
      GVC_PB_X  = PB_X_s,
      GVC_TSImp = TS_imp,
      GVC_TSDom = TS_dom,
      GVC_TS_X  = TS_X_s,
      GVC_X     = GVC_X_s,
      DomX      = DomX_s,
      TradX     = TradX_s,
      X_total   = X_s_tot
    )
  }

  out <- do.call(rbind, output_list)
  rownames(out) <- NULL
  out
}

#' BM_2025 output-based GVC participation indicators
#'
#' For each country s:
#'  - share_GVC_output = GVC_X / X_total
#'  - share_PF_output  = GVC_PF_X / GVC_X
#'  - share_TS_output  = GVC_TS_X / GVC_X
#'  - share_PB_output  = GVC_PB_X / GVC_X
#'  - forward_output   = (GVC_PF_X - GVC_PB_X) / GVC_X
#'
#' @param io bm_io object
#'
#' @return data.frame with exporter-level output-based indicators.
#' @export
bm_2025_output_measures <- function(io) {
  stopifnot(inherits(io, "bm_io"))

  comp <- bm_2025_output_components(io)

  comp$share_GVC_output <- comp$GVC_X / comp$X_total
  comp$share_PF_output  <- comp$GVC_PF_X / comp$GVC_X
  comp$share_TS_output  <- comp$GVC_TS_X / comp$GVC_X
  comp$share_PB_output  <- comp$GVC_PB_X / comp$GVC_X
  comp$forward_output   <- (comp$GVC_PF_X - comp$GVC_PB_X) / comp$GVC_X

  comp
}
