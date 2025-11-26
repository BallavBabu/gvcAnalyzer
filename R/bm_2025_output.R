# R/bm_2025_output.R

# --- Internal Helper (Specific to output calc) -------------------------------

bm_2025_Xexp_list <- function(io) {
  G  <- io$G
  N  <- io$N
  GN <- io$GN
  A  <- io$A

  Xexp_list <- vector("list", G)
  for (r in seq_len(G)) {
    idx_r <- bm_idx_country(io, r)
    A_rr <- A[idx_r, idx_r, drop = FALSE]
    L_rr <- solve(Matrix::Diagonal(N) - A_rr)
    e_rstar <- bm_get_e_star(io, r)
    Xexp_list[[r]] <- as.numeric(L_rr %*% e_rstar)
  }
  Xexp_list
}

#' BM_2025 output-based GVC components by exporter
#' @param io A \code{bm_io} object.
#' @return Data frame with output-based GVC components.
#' @export
bm_2025_output_components <- function(io) {
  stopifnot(inherits(io, "bm_io"))
  G <- io$G; N <- io$N

  L_list    <- bm_L_list(io)
  Xexp_list <- bm_2025_Xexp_list(io)
  output_list <- vector("list", G)

  for (s in seq_len(G)) {
    idx_s <- bm_idx_country(io, s)

    v_s <- io$v[idx_s]
    X_s <- io$X[idx_s]
    L_ss <- L_list[[s]]
    A_ss <- io$A[idx_s, idx_s, drop = FALSE]

    Y_ss    <- io$Y[idx_s, s]
    Y_s_tot <- rowSums(io$Y[idx_s, , drop = FALSE])

    # 1. PURE FORWARD (PF)
    q_total <- numeric(N)
    for (r in seq_len(G)) {
      if (r == s) next
      idx_r  <- bm_idx_country(io, r)
      A_sr   <- io$A[idx_s, idx_r, drop = FALSE]
      Xexp_r <- Xexp_list[[r]]

      inner   <- A_sr %*% Xexp_r
      inner_u <- A_ss %*% (L_ss %*% inner)
      q_total <- q_total + as.numeric(inner + inner_u)
    }
    PF_X_s <- sum(v_s * q_total)

    # 2. PURE BACKWARD (PB)
    fva_tot <- 0
    fva_1b  <- 0

    for (j in seq_len(G)) {
      if (j == s) {
        idx_j <- bm_idx_country(io, j)
        v_j   <- io$v[idx_j]
        L_jj  <- L_list[[j]]

        term_j <- numeric(N)
        for (k in seq_len(G)) {
          if (k == j) next
          idx_k <- bm_idx_country(io, k)
          A_jk  <- io$A[idx_j, idx_k, drop = FALSE]
          B_ks  <- io$B[idx_k, idx_s, drop = FALSE]
          term_j <- term_j + A_jk %*% (B_ks %*% Y_s_tot)
        }
        fva_tot <- fva_tot + as.numeric(t(v_j) %*% L_jj %*% term_j)
      } else {
        idx_j <- bm_idx_country(io, j)
        v_j   <- io$v[idx_j]
        L_jj  <- L_list[[j]]

        term_j <- numeric(N)
        for (k in seq_len(G)) {
          if (k == j) next
          idx_k <- bm_idx_country(io, k)
          A_jk  <- io$A[idx_j, idx_k, drop = FALSE]
          B_ks  <- io$B[idx_k, idx_s, drop = FALSE]
          term_j <- term_j + A_jk %*% (B_ks %*% Y_s_tot)
        }
        fva_tot <- fva_tot + as.numeric(t(v_j) %*% L_jj %*% term_j)

        A_js   <- io$A[idx_j, idx_s, drop = FALSE]
        term_1b <- t(v_j) %*% L_jj %*% A_js %*% L_ss %*% Y_ss
        fva_1b  <- fva_1b + sum(term_1b)
      }
    }
    PB_X_s <- fva_tot - fva_1b

    # 3. TWO-SIDED (TS)
    fva_X_tot <- 0
    fva_X_1b  <- 0

    for (j in seq_len(G)) {
      if (j == s) {
        idx_j <- bm_idx_country(io, j)
        v_j   <- io$v[idx_j]
        L_jj  <- L_list[[j]]

        term_j1 <- numeric(N)
        for (k in seq_len(G)) {
          if (k == j) next
          idx_k <- bm_idx_country(io, k)
          A_jk  <- io$A[idx_j, idx_k, drop = FALSE]
          B_ks  <- io$B[idx_k, idx_s, drop = FALSE]
          term_j1 <- term_j1 + A_jk %*% (B_ks %*% X_s)
        }
        fva_X_tot <- fva_X_tot + as.numeric(t(v_j) %*% L_jj %*% term_j1)
      } else {
        idx_j <- bm_idx_country(io, j)
        v_j   <- io$v[idx_j]
        L_jj  <- L_list[[j]]

        term_j1 <- numeric(N)
        for (k in seq_len(G)) {
          if (k == j) next
          idx_k <- bm_idx_country(io, k)
          A_jk  <- io$A[idx_j, idx_k, drop = FALSE]
          B_ks  <- io$B[idx_k, idx_s, drop = FALSE]
          term_j1 <- term_j1 + A_jk %*% (B_ks %*% X_s)
        }
        fva_X_tot <- fva_X_tot + as.numeric(t(v_j) %*% L_jj %*% term_j1)

        A_js <- io$A[idx_j, idx_s, drop = FALSE]
        term_X_1b <- t(v_j) %*% L_jj %*% A_js %*% L_ss %*% A_ss %*% L_ss %*% Y_ss
        fva_X_1b  <- fva_X_1b + sum(term_X_1b)
      }
    }

    TS_imp <- fva_X_tot - fva_X_1b - PB_X_s
    if (TS_imp < 0) TS_imp <- 0

    term_dom <- A_ss %*% q_total
    TS_dom   <- sum(v_s * as.numeric(L_ss %*% term_dom))

    TS_X_s  <- TS_imp + TS_dom
    GVC_X_s <- PF_X_s + PB_X_s + TS_X_s

    # 4. RESIDUALS
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
#' @param io A \code{bm_io} object.
#' @return Data frame with output-based GVC participation measures.
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
