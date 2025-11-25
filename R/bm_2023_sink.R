# R/bm_2023_sink.R
#
# BM_2023 sink-based bilateral decomposition (/s)
# Importer (sink) perspective for exports from s to r

# Internal: index range of country g
bm_idx_country <- function(io, g) {
  io$.idx_country(g)
}

# -------------------------------------------------------------------
#  BM_2023 sink-based bilateral decomposition
# -------------------------------------------------------------------

#' BM_2023 sink-based bilateral decomposition of exports from s to r
#'
#' Decomposes gross exports e_{sr} (from exporter s to importer r) from the
#' importer (sink) perspective using the BM_2023 /s formulation:
#'
#'   u_N × e_{sr} = DVAsink_sr + DDCsink_sr + FVAsink_sr + FDCsink_sr
#'
#' where the split between \eqn{e^{(s/ \to y^*)}_{sr}} (ultimate shipments)
#' and \eqn{e^{(\to e_s^*)}_{sr}} (re-exports of s's exports) follows
#' Equations (12)–(13) in BM_2023, constructed with the exporter-based
#' slash matrix A^/s (and B^/s).
#'
#' @param io bm_io object as returned by \code{bm_build_io()}
#' @param s  exporter (country index or code, e.g. 1 or "China")
#' @param r  importer (country index or code, e.g. 2 or "India")
#'
#' @return A data frame with one row for the pair (s,r):
#'   \itemize{
#'     \item \code{DVAsink_sr} Domestic VA of s in e_{sr} that is finally
#'           absorbed somewhere in the world as ultimate final demand.
#'     \item \code{DDCsink_sr} Double-counting associated with s in e_{sr}
#'           due to the circulation of s's exports along GVC chains.
#'     \item \code{FVAsink_sr} Foreign VA in e_{sr} from t ≠ s.
#'     \item \code{FDCsink_sr} Foreign double-counting in e_{sr}.
#'     \item \code{EX_sr}      Gross exports e_{sr} (scalar).
#'   }
#' @export
bm_2023_bilateral_sink <- function(io, s, r) {
  stopifnot(inherits(io, "bm_io"))
  s <- bm_country_id(io, s)
  r <- bm_country_id(io, r)
  if (s == r) stop("s and r must be different")

  G    <- io$G
  N    <- io$N
  GN   <- io$GN
  vvec <- io$v
  A    <- io$A
  B    <- io$B
  Y    <- io$Y

  idx_s <- bm_idx_country(io, s)
  idx_r <- bm_idx_country(io, r)

  v_s <- vvec[idx_s]

  # ------------------------------------------------------------------
  # Basic blocks: B_ss, A_sr, A_rr, L_rr
  # ------------------------------------------------------------------
  B_ss <- B[idx_s, idx_s, drop = FALSE]

  A_sr <- A[idx_s, idx_r, drop = FALSE]
  A_rr <- A[idx_r, idx_r, drop = FALSE]
  I_N  <- Matrix::Diagonal(n = N)
  L_rr <- solve(I_N - A_rr)   # domestic Leontief for r

  # Final-demand pieces
  y_sr <- Y[idx_s, r]                          # s -> r final demand
  y_rr <- Y[idx_r, r]                          # r -> r final demand

  # Σ_{j≠r} y_rj (all destinations of r except r itself)
  sum_y_r_other <- rep(0, N)
  for (j in seq_len(G)) {
    if (j == r) next
    sum_y_r_other <- sum_y_r_other + Y[idx_r, j]
  }

  # ------------------------------------------------------------------
  # Split of e_{sr} into:
  #   e_sy_sr = e^{(s/ -> y^*)}_{sr}   (ultimate-shipments component)
  #   e_es_sr = e^{(-> e_s^*)}_{sr}    (re-export / re-routing component)
  #   following BM_2023 eqs. (12)–(13) using B^/s
  # ------------------------------------------------------------------

  # B^/s: exporter-based slash Leontief (intermediate exports of s removed)
  Bs <- bm_Bslash_s(io, s)

  # Precompute Y totals by origin k: y_k* = sum_l y_kl
  Y_total_list <- vector("list", G)
  for (k in seq_len(G)) {
    Y_total_list[[k]] <- rowSums(Y[bm_idx_country(io, k), , drop = FALSE])
  }

  # Country-total exports from s: e_s* = sum_{z≠s} e_sz
  e_s_star <- Reduce(
    `+`,
    lapply(setdiff(seq_len(G), s), function(rr) bm_e_sr(io, s, rr))
  )

  # x(/s -> y*)_j and x(-> e_s*)_j, as in BM eq. (12)–(13)
  x_sy_for_j <- function(j) {
    idx_j <- bm_idx_country(io, j)

    # x(/s -> y*)_j = Σ_{k≠s} B^/s_jk y_k* + B^/s_js y_ss
    x_j <- rep(0, N)

    for (k in seq_len(G)) {
      if (k == s) next
      idx_k <- bm_idx_country(io, k)
      B_jk  <- Bs[idx_j, idx_k, drop = FALSE]
      y_k_total <- Y_total_list[[k]]
      x_j <- x_j + B_jk %*% y_k_total
    }

    # Add B^/s_js y_ss (domestic final demand of s)
    idx_s_local <- bm_idx_country(io, s)
    B_js  <- Bs[idx_j, idx_s_local, drop = FALSE]
    y_ss  <- Y[idx_s_local, s]
    x_j   <- x_j + B_js %*% y_ss

    as.numeric(x_j)
  }

  x_es_for_j <- function(j) {
    idx_j <- bm_idx_country(io, j)
    idx_s_local <- bm_idx_country(io, s)
    B_js  <- Bs[idx_j, idx_s_local, drop = FALSE]
    as.numeric(B_js %*% e_s_star)
  }

  # Base terms in BM eq. (12) that belong entirely to "ultimate shipments"
  base_ult <- y_sr +
    A_sr %*% (L_rr %*% y_rr) +
    A_sr %*% (L_rr %*% sum_y_r_other)

  # Decompose the last term in eq. (12) into /s->y* and ->e_s*
  sum_A_rj_x_sy <- rep(0, N)
  sum_A_rj_x_es <- rep(0, N)

  for (j in seq_len(G)) {
    if (j == r) next
    idx_j <- bm_idx_country(io, j)
    A_rj  <- A[idx_r, idx_j, drop = FALSE]

    x_sy_j <- x_sy_for_j(j)
    x_es_j <- x_es_for_j(j)

    sum_A_rj_x_sy <- sum_A_rj_x_sy + A_rj %*% x_sy_j
    sum_A_rj_x_es <- sum_A_rj_x_es + A_rj %*% x_es_j
  }

  extra_ult <- A_sr %*% (L_rr %*% sum_A_rj_x_sy)
  extra_es  <- A_sr %*% (L_rr %*% sum_A_rj_x_es)

  # Components of e_sr
  e_sy_sr <- as.numeric(base_ult + extra_ult)  # N-vector
  e_es_sr <- as.numeric(extra_es)              # N-vector

  # Consistency: total e_sr
  e_sr_vec <- bm_e_sr(io, s, r)
  EX_sr    <- sum(e_sr_vec)

  # ------------------------------------------------------------------
  # Map to DVAsink, DDCsink, FVAsink, FDCsink (BM 4-way decomposition)
  # ------------------------------------------------------------------

  # DVAsink_sr = v_s × B_ss × e_sy_sr
  vec_DVAsink <- B_ss %*% e_sy_sr
  DVAsink     <- sum(v_s * as.numeric(vec_DVAsink))

  # DDCsink_sr = v_s × B_ss × e_es_sr
  vec_DDCsink <- B_ss %*% e_es_sr
  DDCsink     <- sum(v_s * as.numeric(vec_DDCsink))

  # FVAsink_sr = Σ_{t≠s} v_t × B_ts × e_sy_sr
  # FDCsink_sr = Σ_{t≠s} v_t × B_ts × e_es_sr
  FVAsink <- 0
  FDCsink <- 0
  for (t in seq_len(G)) {
    if (t == s) next
    idx_t <- bm_idx_country(io, t)
    v_t   <- vvec[idx_t]
    B_ts  <- B[idx_t, idx_s, drop = FALSE]

    vec_FVA_t <- B_ts %*% e_sy_sr
    FVAsink   <- FVAsink + sum(v_t * as.numeric(vec_FVA_t))

    vec_FDC_t <- B_ts %*% e_es_sr
    FDCsink   <- FDCsink + sum(v_t * as.numeric(vec_FDC_t))
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
