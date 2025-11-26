# R/bm_2023_sink.R

#' BM_2023 sink-based bilateral decomposition of exports from s to r
#'
#' @param io A \code{bm_io} object.
#' @param s Exporter country (name or index).
#' @param r Importer country (name or index).
#' @return A data frame with the sink-based value-added decomposition.
#' @export
bm_2023_bilateral_sink <- function(io, s, r) {
  stopifnot(inherits(io, "bm_io"))
  s <- bm_country_id(io, s)
  r <- bm_country_id(io, r)
  if (s == r) stop("s and r must be different")

  G    <- io$G
  N    <- io$N
  vvec <- io$v
  A    <- io$A
  B    <- io$B
  Y    <- io$Y

  idx_s <- bm_idx_country(io, s)
  idx_r <- bm_idx_country(io, r)
  v_s <- vvec[idx_s]

  B_ss <- bm_block(io, B, s, s)
  A_sr <- bm_block(io, A, s, r)
  A_rr <- bm_block(io, A, r, r)

  I_N  <- Matrix::Diagonal(n = N)
  L_rr <- solve(I_N - A_rr)

  y_sr <- Y[idx_s, r]
  y_rr <- Y[idx_r, r]

  sum_y_r_other <- rep(0, N)
  for (j in seq_len(G)) {
    if (j == r) next
    sum_y_r_other <- sum_y_r_other + Y[idx_r, j]
  }

  Bs <- bm_Bslash_s(io, s)

  Y_total_list <- vector("list", G)
  for (k in seq_len(G)) {
    Y_total_list[[k]] <- rowSums(Y[bm_idx_country(io, k), , drop = FALSE])
  }

  e_s_star <- bm_get_e_star(io, s)

  x_sy_for_j <- function(j) {
    idx_j <- bm_idx_country(io, j)
    x_j <- rep(0, N)

    for (k in seq_len(G)) {
      if (k == s) next
      B_jk  <- bm_block(io, Bs, j, k)
      y_k_total <- Y_total_list[[k]]
      x_j <- x_j + B_jk %*% y_k_total
    }

    idx_s_local <- bm_idx_country(io, s)
    B_js  <- bm_block(io, Bs, j, s)
    y_ss  <- Y[idx_s_local, s]
    x_j   <- x_j + B_js %*% y_ss

    as.numeric(x_j)
  }

  x_es_for_j <- function(j) {
    B_js  <- bm_block(io, Bs, j, s)
    as.numeric(B_js %*% e_s_star)
  }

  base_ult <- y_sr +
    A_sr %*% (L_rr %*% y_rr) +
    A_sr %*% (L_rr %*% sum_y_r_other)

  sum_A_rj_x_sy <- rep(0, N)
  sum_A_rj_x_es <- rep(0, N)

  for (j in seq_len(G)) {
    if (j == r) next
    A_rj  <- bm_block(io, A, r, j)

    x_sy_j <- x_sy_for_j(j)
    x_es_j <- x_es_for_j(j)

    sum_A_rj_x_sy <- sum_A_rj_x_sy + A_rj %*% x_sy_j
    sum_A_rj_x_es <- sum_A_rj_x_es + A_rj %*% x_es_j
  }

  extra_ult <- A_sr %*% (L_rr %*% sum_A_rj_x_sy)
  extra_es  <- A_sr %*% (L_rr %*% sum_A_rj_x_es)

  e_sy_sr <- as.numeric(base_ult + extra_ult)
  e_es_sr <- as.numeric(extra_es)

  e_sr_vec <- bm_get_e_sr(io, s, r)
  EX_sr    <- sum(e_sr_vec)

  DVAsink     <- sum(v_s * as.numeric(B_ss %*% e_sy_sr))
  DDCsink     <- sum(v_s * as.numeric(B_ss %*% e_es_sr))

  FVAsink <- 0
  FDCsink <- 0
  for (t in seq_len(G)) {
    if (t == s) next
    v_t   <- vvec[bm_idx_country(io, t)]
    B_ts  <- bm_block(io, B, t, s)

    FVAsink   <- FVAsink + sum(v_t * as.numeric(B_ts %*% e_sy_sr))
    FDCsink   <- FDCsink + sum(v_t * as.numeric(B_ts %*% e_es_sr))
  }

  data.frame(
    exporter    = io$countries[s],
    importer    = io$countries[r],
    DVAsink_sr  = DVAsink,
    DDCsink_sr  = DDCsink,
    FVAsink_sr  = FVAsink,
    FDCsink_sr  = FDCsink,
    EX_sr       = EX_sr
  )
}
