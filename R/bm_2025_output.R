# R/bm_2025_output.R
#
# BM_2025 output-based GVC measures (PF/TS/PB) and indicators.
# FIXED:
#   (i) Xexp_r now uses domestic Leontief L_rr (not B_rr),
#   (ii) PF = v_s * q_total (no extra L_ss), consistent with BMT (2021).

# --- Local Helpers -----------------------------------------------------------
bm25_idx_country <- function(io, g) {
  if (is.character(g)) g <- match(g, io$countries)
  ((g - 1L) * io$N + 1L):(g * io$N)
}

bm25_L_list <- function(io) {
  G <- io$G; N <- io$N; A <- io$A
  L_list <- vector("list", G)
  for (g in seq_len(G)) {
    idx_g <- bm25_idx_country(io, g)
    A_gg  <- A[idx_g, idx_g, drop = FALSE]
    L_list[[g]] <- solve(Matrix::Diagonal(N) - A_gg)
  }
  L_list
}

bm_2025_Xexp_list <- function(io) {
  # NOTE: Xexp_r = L_rr * sum_{k != r} E_rk  (BMT, 2021).
  G  <- io$G
  N  <- io$N
  GN <- io$GN
  A  <- io$A

  Xexp_list <- vector("list", G)
  for (r in seq_len(G)) {
    idx_r <- bm25_idx_country(io, r)

    # Domestic Leontief for country r
    A_rr <- A[idx_r, idx_r, drop = FALSE]
    L_rr <- solve(Matrix::Diagonal(N) - A_rr)

    # Exports from r to all k != r (intermediate + final)
    int_exp <- rowSums(io$Z[idx_r, setdiff(seq_len(GN), idx_r), drop = FALSE])
    fin_exp <- rowSums(io$Y[idx_r, setdiff(seq_len(G),  r),   drop = FALSE])
    e_rstar <- int_exp + fin_exp

    Xexp_list[[r]] <- as.numeric(L_rr %*% e_rstar)
  }
  Xexp_list
}

#' BM_2025 output-based GVC components by exporter
#' @export
bm_2025_output_components <- function(io) {
  stopifnot(inherits(io, "bm_io"))
  G <- io$G; N <- io$N

  L_list    <- bm25_L_list(io)
  Xexp_list <- bm_2025_Xexp_list(io)
  output_list <- vector("list", G)

  for (s in seq_len(G)) {
    idx_s <- bm25_idx_country(io, s)

    v_s <- io$v[idx_s]              # value-added coefficients (length N)
    X_s <- io$X[idx_s]              # gross output by sector of s
    L_ss <- L_list[[s]]             # domestic Leontief
    A_ss <- io$A[idx_s, idx_s, drop = FALSE]

    Y_ss    <- io$Y[idx_s, s]                     # final demand of s in s
    Y_s_tot <- rowSums(io$Y[idx_s, , drop = FALSE])  # total final demand of s

    ## 1. PURE FORWARD (PF) -----------------------------------------------
    # q_total = Σ_{r != s} [ A_sr Xexp_r + A_ss L_ss A_sr Xexp_r ]
    # PF_X_s  = v_s * q_total  (GVC value-added traced in output)
    q_total <- numeric(N)
    for (r in seq_len(G)) {
      if (r == s) next
      idx_r  <- bm25_idx_country(io, r)
      A_sr   <- io$A[idx_s, idx_r, drop = FALSE]
      Xexp_r <- Xexp_list[[r]]

      inner   <- A_sr %*% Xexp_r                # A_sr Xexp_r
      inner_u <- A_ss %*% (L_ss %*% inner)      # A_ss L_ss A_sr Xexp_r
      q_total <- q_total + as.numeric(inner + inner_u)
    }
    # CORRECT: PF = value added of q_total (no extra L_ss)
    PF_X_s <- sum(v_s * q_total)

    ## 2. PURE BACKWARD (PB) ----------------------------------------------
    # Follows BMT (2021) "GVC-output: Pure Backward" logic.
    fva_tot <- 0
    fva_1b  <- 0

    for (j in seq_len(G)) {
      if (j == s) {
        # j = s contributes only to the first term (k != j restriction handled below)
        idx_j <- bm25_idx_country(io, j)
        v_j   <- io$v[idx_j]
        L_jj  <- L_list[[j]]

        term_j <- numeric(N)
        for (k in seq_len(G)) {
          if (k == j) next
          idx_k <- bm25_idx_country(io, k)
          A_jk  <- io$A[idx_j, idx_k, drop = FALSE]
          B_ks  <- io$B[idx_k, idx_s, drop = FALSE]
          term_j <- term_j + A_jk %*% (B_ks %*% Y_s_tot)
        }
        fva_tot <- fva_tot + as.numeric(t(v_j) %*% L_jj %*% term_j)
      } else {
        # j != s: contributes to both terms
        idx_j <- bm25_idx_country(io, j)
        v_j   <- io$v[idx_j]
        L_jj  <- L_list[[j]]

        # First term: Σ_{k != j} A_jk B_ks Y_s_tot
        term_j <- numeric(N)
        for (k in seq_len(G)) {
          if (k == j) next
          idx_k <- bm25_idx_country(io, k)
          A_jk  <- io$A[idx_j, idx_k, drop = FALSE]
          B_ks  <- io$B[idx_k, idx_s, drop = FALSE]
          term_j <- term_j + A_jk %*% (B_ks %*% Y_s_tot)
        }
        fva_tot <- fva_tot + as.numeric(t(v_j) %*% L_jj %*% term_j)

        # Second term: Σ_{j != s} V_j L_jj A_js L_ss Y_ss
        A_js   <- io$A[idx_j, idx_s, drop = FALSE]
        term_1b <- t(v_j) %*% L_jj %*% A_js %*% L_ss %*% Y_ss
        fva_1b  <- fva_1b + sum(term_1b)
      }
    }
    PB_X_s <- fva_tot - fva_1b

    ## 3. TWO-SIDED (TS) --------------------------------------------------
    # Imported inputs re-exported as intermediates:
    # GVCTwoSidedImpInp_s =
    #   Σ_j V_j L_jj Σ_{k != j} A_jk B_ks X_s
    #   - Σ_{j != s} V_j L_jj A_js L_ss A_ss L_ss Y_ss
    #   - GVCPureBackX_s
    fva_X_tot <- 0
    fva_X_1b  <- 0

    for (j in seq_len(G)) {
      if (j == s) {
        idx_j <- bm25_idx_country(io, j)
        v_j   <- io$v[idx_j]
        L_jj  <- L_list[[j]]

        term_j1 <- numeric(N)
        for (k in seq_len(G)) {
          if (k == j) next
          idx_k <- bm25_idx_country(io, k)
          A_jk  <- io$A[idx_j, idx_k, drop = FALSE]
          B_ks  <- io$B[idx_k, idx_s, drop = FALSE]
          term_j1 <- term_j1 + A_jk %*% (B_ks %*% X_s)
        }
        fva_X_tot <- fva_X_tot + as.numeric(t(v_j) %*% L_jj %*% term_j1)
      } else {
        idx_j <- bm25_idx_country(io, j)
        v_j   <- io$v[idx_j]
        L_jj  <- L_list[[j]]

        # First term
        term_j1 <- numeric(N)
        for (k in seq_len(G)) {
          if (k == j) next
          idx_k <- bm25_idx_country(io, k)
          A_jk  <- io$A[idx_j, idx_k, drop = FALSE]
          B_ks  <- io$B[idx_k, idx_s, drop = FALSE]
          term_j1 <- term_j1 + A_jk %*% (B_ks %*% X_s)
        }
        fva_X_tot <- fva_X_tot + as.numeric(t(v_j) %*% L_jj %*% term_j1)

        # Second term for two-sided imported inputs:
        # Σ_{j != s} V_j L_jj A_js L_ss A_ss L_ss Y_ss
        A_js <- io$A[idx_j, idx_s, drop = FALSE]
        term_X_1b <- t(v_j) %*% L_jj %*% A_js %*% L_ss %*% A_ss %*% L_ss %*% Y_ss
        fva_X_1b  <- fva_X_1b + sum(term_X_1b)
      }
    }

    TS_imp <- fva_X_tot - fva_X_1b - PB_X_s
    if (TS_imp < 0) TS_imp <- 0

    # Domestic inputs sold along GVC chains:
    # GVCTwoSidedDomInp_s = V_s L_ss A_ss Σ_{r != s} (A_sr Xexp_r + A_ss L_ss A_sr Xexp_r)
    term_dom <- A_ss %*% q_total
    TS_dom   <- sum(v_s * as.numeric(L_ss %*% term_dom))

    TS_X_s  <- TS_imp + TS_dom
    GVC_X_s <- PF_X_s + PB_X_s + TS_X_s

    ## 4. RESIDUALS --------------------------------------------------------
    DomX_s  <- sum(v_s * as.numeric(L_ss %*% Y_ss))
    TradX_s <- sum(X_s) - DomX_s - GVC_X_s

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
      X_total   = sum(X_s)
    )
  }

  do.call(rbind, output_list)
}

#' BM_2025 output-based GVC participation indicators
#' @export
bm_2025_output_measures <- function(io) {
  stopifnot(inherits(io, "bm_io"))
  comp <- bm_2025_output_components(io)

  comp$share_GVC_output <- comp$GVC_X   / comp$X_total
  comp$share_PF_output  <- ifelse(comp$GVC_X > 0, comp$GVC_PF_X / comp$GVC_X, 0)
  comp$share_TS_output  <- ifelse(comp$GVC_X > 0, comp$GVC_TS_X / comp$GVC_X, 0)
  comp$share_PB_output  <- ifelse(comp$GVC_X > 0, comp$GVC_PB_X / comp$GVC_X, 0)
  comp$forward_output   <- ifelse(comp$GVC_X > 0,
                                  (comp$GVC_PF_X - comp$GVC_PB_X) / comp$GVC_X, 0)
  comp
}
