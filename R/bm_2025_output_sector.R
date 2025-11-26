# R/bm_2025_output_sector.R

#' BM 2025 output components by country and sector
#' @param io A \code{bm_io} object.
#' @return Data frame with sectoral output components.
#' @export
bm_2025_output_components_sector <- function(io) {
  stopifnot(inherits(io, "bm_io"))

  G <- io$G; N <- io$N; GN <- io$GN
  A <- io$A; B <- io$B; Y <- io$Y; X <- io$X; v <- io$v

  L_list <- bm_L_list(io)

  Xexp_list <- vector("list", G)
  for (r in seq_len(G)) {
    e_r_star <- bm_get_e_star(io, r)
    L_rr <- L_list[[r]]
    Xexp_list[[r]] <- as.numeric(L_rr %*% e_r_star)
  }

  res_list <- vector("list", G)

  for (s in seq_len(G)) {
    idx_s <- bm_idx_country(io, s)
    v_s_vec <- v[idx_s]
    X_s_vec <- X[idx_s]
    L_ss    <- L_list[[s]]
    A_ss    <- A[idx_s, idx_s, drop = FALSE]
    Y_ss    <- Y[idx_s, s]
    Y_s_tot <- rowSums(Y[idx_s, , drop = FALSE])

    # 1. PURE FORWARD (PF)
    q_total <- numeric(N)
    for (r in seq_len(G)) {
      if (r == s) next
      idx_r  <- bm_idx_country(io, r)
      A_sr   <- A[idx_s, idx_r, drop = FALSE]
      Xexp_r <- Xexp_list[[r]]

      inner   <- A_sr %*% Xexp_r
      inner_u <- A_ss %*% (L_ss %*% inner)
      q_total <- q_total + as.numeric(inner + inner_u)
    }
    PF_vec <- v_s_vec * q_total

    # 2. PURE BACKWARD (PB)
    fva_intensity_total <- numeric(N)
    fva_intensity_1border <- numeric(N)

    for (j in seq_len(G)) {
      if (j == s) next
      idx_j <- bm_idx_country(io, j)
      v_j   <- v[idx_j]
      B_js  <- bm_block(io, B, j, s)
      fva_intensity_total <- fva_intensity_total + as.numeric(t(v_j) %*% B_js)

      L_jj <- L_list[[j]]
      A_js <- bm_block(io, A, j, s)

      term_1b <- as.numeric(t(v_j) %*% L_jj %*% A_js %*% L_ss)
      fva_intensity_1border <- fva_intensity_1border + term_1b
    }
    PB_vec <- (fva_intensity_total * Y_s_tot) - (fva_intensity_1border * Y_ss)

    # 3. Two-Sided (TS)
    TS_Imp_vec <- (fva_intensity_total * X_s_vec) - (fva_intensity_1border * Y_ss) - PB_vec
    TS_Imp_vec[TS_Imp_vec < 0] <- 0

    term_dom <- A_ss %*% q_total
    TS_Dom_vec <- v_s_vec * as.numeric(L_ss %*% term_dom)

    TS_vec <- TS_Imp_vec + TS_Dom_vec
    GVC_vec <- PF_vec + PB_vec + TS_vec

    # 4. Residuals
    DomX_vec <- v_s_vec * as.numeric(L_ss %*% Y_ss)
    TradX_vec <- X_s_vec - DomX_vec - GVC_vec

    res_list[[s]] <- data.frame(
      country      = rep(io$countries[s], N),
      sector       = io$sectors,
      X_i          = as.numeric(X_s_vec),
      DomX_i       = as.numeric(DomX_vec),
      TradX_i      = as.numeric(TradX_vec),
      GVC_PF_Xi    = as.numeric(PF_vec),
      GVC_PB_Xi    = as.numeric(PB_vec),
      GVC_TSImp_i  = as.numeric(TS_Imp_vec),
      GVC_TSDom_i  = as.numeric(TS_Dom_vec),
      GVC_TS_Xi    = as.numeric(TS_vec),
      GVC_Xi       = as.numeric(GVC_vec)
    )
  }
  do.call(rbind, res_list)
}

#' BM 2025 output participation measures by country and sector
#' @param io A \code{bm_io} object.
#' @return Data frame with sectoral GVC measures.
#' @export
bm_2025_output_measures_sector <- function(io) {
  stopifnot(inherits(io, "bm_io"))
  df <- bm_2025_output_components_sector(io)

  df$share_GVC_output_i <- ifelse(df$X_i > 0, df$GVC_Xi / df$X_i, NA_real_)
  df$share_PF_output_i  <- ifelse(df$GVC_Xi > 0, df$GVC_PF_Xi / df$GVC_Xi, NA_real_)
  df$share_TS_output_i  <- ifelse(df$GVC_Xi > 0, df$GVC_TS_Xi / df$GVC_Xi, NA_real_)
  df$share_PB_output_i  <- ifelse(df$GVC_Xi > 0, df$GVC_PB_Xi / df$GVC_Xi, NA_real_)
  df$forward_output_i   <- ifelse(df$GVC_Xi > 0, (df$GVC_PF_Xi - df$GVC_PB_Xi) / df$GVC_Xi, NA_real_)
  df
}
