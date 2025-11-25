# bm_gvc_measures.R
# Summary GVC participation indicators (BM_2025 Sections 3.3 and 4)

# This file assumes:
# - bm_tripartite_trade_all(io)   ## from bm_tripartite_trade.R
# - bm_gvc_output_all(io)         ## from bm_output_gvc.R

# ------------------------------------------------------------
# 4.1 Trade-based GVC participation indicators
# ------------------------------------------------------------

#' Trade-based GVC decomposition by exporter
#'
#' Aggregates the bilateral tripartite GVC trade decomposition
#' over all importing partners for each exporting country s.
#'
#' For each exporter, it returns:
#' \itemize{
#'   \item E_sr     : total gross exports (sum over r)
#'   \item GVC_sr   : total GVC-related trade
#'   \item GVC_PF   : pure-forward GVC-related trade
#'   \item GVC_TS   : two-sided GVC-related trade
#'   \item GVC_PB   : pure-backward GVC-related trade
#' }
#'
#' @param io A \code{bm_io} object.
#'
#' @return A \code{data.frame} with one row per exporter.
#' @export
bm_gvc_trade_by_exporter <- function(io) {
  trip_df <- bm_tripartite_trade_all(io)
  # aggregate by exporter (country name)
  agg <- aggregate(
    cbind(E_sr, GVC_sr, GVC_PF, GVC_TS, GVC_PB) ~ exporter,
    data = trip_df,
    FUN = sum
  )
  # keep a clean country column name
  names(agg)[names(agg) == "exporter"] <- "country"
  agg
}

#' Trade-based GVC participation indicators
#'
#' Builds on \code{bm_gvc_trade_by_exporter()} and computes,
#' for each exporter:
#' \itemize{
#'   \item share_GVC_trade = GVC_sr / E_sr
#'   \item share_PF_trade  = GVC_PF / GVC_sr
#'   \item share_TS_trade  = GVC_TS / GVC_sr
#'   \item share_PB_trade  = GVC_PB / GVC_sr
#'   \item forward_trade   = (GVC_PF - GVC_PB) / GVC_sr
#' }
#'
#' These correspond to the paper's trade-based participation
#' and forwardness indicators.
#'
#' @param io A \code{bm_io} object.
#'
#' @return A \code{data.frame} with one row per country.
#' @export
bm_gvc_trade_measures <- function(io) {
  df <- bm_gvc_trade_by_exporter(io)

  within(df, {
    share_GVC_trade <- ifelse(E_sr > 0, GVC_sr / E_sr, NA_real_)
    share_PF_trade  <- ifelse(GVC_sr > 0, GVC_PF / GVC_sr, NA_real_)
    share_TS_trade  <- ifelse(GVC_sr > 0, GVC_TS / GVC_sr, NA_real_)
    share_PB_trade  <- ifelse(GVC_sr > 0, GVC_PB / GVC_sr, NA_real_)
    forward_trade   <- ifelse(GVC_sr > 0, (GVC_PF - GVC_PB) / GVC_sr, NA_real_)
  })
}

# ------------------------------------------------------------
# 4.2 Output-based GVC participation indicators
# ------------------------------------------------------------

#' Output-based GVC participation indicators
#'
#' Builds on \code{bm_gvc_output_all()} and computes, for each country:
#' \itemize{
#'   \item share_GVC_output = GVC_X / X_total
#'   \item share_PF_output  = GVC_PF_X / GVC_X
#'   \item share_TS_output  = GVC_TS_X / GVC_X
#'   \item share_PB_output  = GVC_PB_X / GVC_X
#'   \item forward_output   = (GVC_PF_X - GVC_PB_X) / GVC_X
#' }
#'
#' These mirror the trade-based indicators on the output side.
#'
#' @param io A \code{bm_io} object.
#'
#' @return A \code{data.frame} with one row per country.
#' @export
bm_gvc_output_measures <- function(io) {
  df <- bm_gvc_output_all(io)

  within(df, {
    share_GVC_output <- ifelse(X_total > 0, GVC_X / X_total, NA_real_)
    share_PF_output  <- ifelse(GVC_X > 0, GVC_PF_X / GVC_X, NA_real_)
    share_TS_output  <- ifelse(GVC_X > 0, GVC_TS_X / GVC_X, NA_real_)
    share_PB_output  <- ifelse(GVC_X > 0, GVC_PB_X / GVC_X, NA_real_)
    forward_output   <- ifelse(GVC_X > 0, (GVC_PF_X - GVC_PB_X) / GVC_X, NA_real_)
  })
}

# ------------------------------------------------------------
# 4.3 Convenience wrapper: trade + output together
# ------------------------------------------------------------

#' Combined trade- and output-based GVC measures
#'
#' Returns a list with:
#' \itemize{
#'   \item trade  : trade-based decomposition + indicators
#'   \item output : output-based decomposition + indicators
#' }
#'
#' @param io A \code{bm_io} object.
#'
#' @return A named list with components \code{trade} and \code{output}.
#' @export
bm_gvc_measures_all <- function(io) {
  list(
    trade  = bm_gvc_trade_measures(io),
    output = bm_gvc_output_measures(io)
  )
}
