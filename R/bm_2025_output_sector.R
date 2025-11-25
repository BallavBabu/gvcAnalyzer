# R/bm_2025_output_sector.R
#
# BM 2025 output indicators at country × sector level.
# This extends the country-level BM 2025 output decomposition
# by allocating each country's BM 2025 components to sectors
# in proportion to sectoral gross output X_{si}.
#
# This preserves country totals from bm_2025_output_components(io).

#' BM 2025 output components by country and sector
#'
#' This function extends \code{bm_2025_output_components()} to the
#' country–sector level. For each country \eqn{s} and sector \eqn{i},
#' it allocates the country-level BM 2025 output components
#' (DomX, TradX, GVC_PF_X, GVC_PB_X, GVC_TSImp, GVC_TSDom, GVC_TS_X, GVC_X)
#' in proportion to the sector's share in total gross output:
#'
#' \deqn{
#'   \theta_{si} = X_{si} / \sum_{j} X_{sj},
#' }
#'
#' and
#'
#' \deqn{
#'   \mathrm{DomX}_{si}    = \theta_{si} \cdot \mathrm{DomX}_s, \quad
#'   \mathrm{TradX}_{si}   = \theta_{si} \cdot \mathrm{TradX}_s, \quad
#'   \mathrm{GVC\_PF\_X}_{si}  = \theta_{si} \cdot \mathrm{GVC\_PF\_X}_s,
#' }
#'
#' etc. This ensures that for each country \eqn{s}:
#'
#' \deqn{
#'   \sum_i \mathrm{DomX}_{si}    = \mathrm{DomX}_s, \quad
#'   \sum_i \mathrm{TradX}_{si}   = \mathrm{TradX}_s, \quad
#'   \sum_i \mathrm{GVC\_PF\_X}_{si} = \mathrm{GVC\\_PF\_X}_s,
#' }
#'
#' and similarly for the other BM 2025 components.
#'
#' @param io A \code{bm_io} object created by \code{bm_build_io()}.
#'
#' @return A data frame with one row per country–sector, containing:
#'   \itemize{
#'     \item \code{country}: country name
#'     \item \code{sector}: sector name
#'     \item \code{X_i}: sectoral gross output \eqn{X_{si}}
#'     \item \code{DomX_i}, \code{TradX_i}
#'     \item \code{GVC_PF_Xi}, \code{GVC_PB_Xi}
#'     \item \code{GVC_TSImp_i}, \code{GVC_TSDom_i}
#'     \item \code{GVC_TS_Xi}, \code{GVC_Xi}
#'   }
#' @export
bm_2025_output_components_sector <- function(io) {
  stopifnot(inherits(io, "bm_io"))

  # country-level BM 2025 output components
  comp_country <- bm_2025_output_components(io)

  res_list <- vector("list", io$G * io$N)
  idx <- 1L

  for (g in seq_len(io$G)) {
    c_name <- io$countries[g]

    row_c <- comp_country[comp_country$country == c_name, , drop = FALSE]
    if (nrow(row_c) != 1L) {
      stop("bm_2025_output_components: country row not found for ", c_name)
    }

    # Sectoral outputs X_{si} for this country
    idx_g <- bm_idx_country(io, g)
    X_g   <- io$X[idx_g]
    X_tot <- sum(X_g)

    if (X_tot == 0) {
      # if a country has zero total output, use equal weights
      share <- rep(1 / io$N, io$N)
    } else {
      share <- X_g / X_tot
    }

    for (i in seq_len(io$N)) {
      res_list[[idx]] <- data.frame(
        country      = c_name,
        sector       = io$sectors[i],
        X_i          = X_g[i],
        DomX_i       = row_c$DomX       * share[i],
        TradX_i      = row_c$TradX      * share[i],
        GVC_PF_Xi    = row_c$GVC_PF_X   * share[i],
        GVC_PB_Xi    = row_c$GVC_PB_X   * share[i],
        GVC_TSImp_i  = row_c$GVC_TSImp  * share[i],
        GVC_TSDom_i  = row_c$GVC_TSDom  * share[i],
        GVC_TS_Xi    = row_c$GVC_TS_X   * share[i],
        GVC_Xi       = row_c$GVC_X      * share[i]
      )
      idx <- idx + 1L
    }
  }

  out <- do.call(rbind, res_list)
  rownames(out) <- NULL
  out
}

#' BM 2025 output participation measures by country and sector
#'
#' Based on \code{bm_2025_output_components_sector()}, this function computes
#' sector-level BM 2025 participation measures:
#'
#' \deqn{
#'   \mathrm{share\_GVC\_output}_{si} = \mathrm{GVC\_X}_{si} / X_{si},
#' }
#'
#' \deqn{
#'   \mathrm{share\_PF\_output}_{si} = \mathrm{GVC\_PF\_X}_{si} / \mathrm{GVC\_X}_{si},
#' }
#'
#' \deqn{
#'   \mathrm{share\_TS\_output}_{si} = \mathrm{GVC\_TS\_X}_{si} / \mathrm{GVC\_X}_{si},
#' }
#'
#' \deqn{
#'   \mathrm{share\_PB\_output}_{si} = \mathrm{GVC\_PB\_X}_{si} / \mathrm{GVC\_X}_{si},
#' }
#'
#' and a sector-level forwardness index:
#'
#' \deqn{
#'   \mathrm{forward\_output}_{si}
#'   = \left(\mathrm{GVC\_PF\_X}_{si} - \mathrm{GVC\_PB\_X}_{si}\right)
#'     / \mathrm{GVC\_X}_{si}.
#' }
#'
#' Sectors with zero output or zero GVC-related output receive \code{NA}
#' for the corresponding ratios.
#'
#' @param io A \code{bm_io} object created by \code{bm_build_io()}.
#'
#' @return A data frame with the columns from
#'   \code{bm_2025_output_components_sector()} plus:
#'   \itemize{
#'     \item \code{share_GVC_output_i}
#'     \item \code{share_PF_output_i}
#'     \item \code{share_TS_output_i}
#'     \item \code{share_PB_output_i}
#'     \item \code{forward_output_i}
#'   }
#' @export
bm_2025_output_measures_sector <- function(io) {
  stopifnot(inherits(io, "bm_io"))

  df <- bm_2025_output_components_sector(io)

  df$share_GVC_output_i <- ifelse(df$X_i > 0,
                                  df$GVC_Xi / df$X_i,
                                  NA_real_)
  df$share_PF_output_i  <- ifelse(df$GVC_Xi > 0,
                                  df$GVC_PF_Xi / df$GVC_Xi,
                                  NA_real_)
  df$share_TS_output_i  <- ifelse(df$GVC_Xi > 0,
                                  df$GVC_TS_Xi / df$GVC_Xi,
                                  NA_real_)
  df$share_PB_output_i  <- ifelse(df$GVC_Xi > 0,
                                  df$GVC_PB_Xi / df$GVC_Xi,
                                  NA_real_)
  df$forward_output_i   <- ifelse(df$GVC_Xi > 0,
                                  (df$GVC_PF_Xi - df$GVC_PB_Xi) / df$GVC_Xi,
                                  NA_real_)

  df
}
