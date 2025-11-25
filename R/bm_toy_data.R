# R/bm_toy_data.R

#' Toy 4-country, 3-sector IO table for bmGVC
#'
#' A small multi-country input–output data set used in bmGVC examples
#' and vignettes. It contains four countries (China, India, Japan, ROW)
#' and three sectors (Primary, Manufacturing, Service).
#'
#' The data are stored in six objects:
#' \itemize{
#'   \item \code{bm_toy_Z}: 12 x 12 intermediate demand matrix
#'   \item \code{bm_toy_Y}: 12 x 4 final demand matrix
#'   \item \code{bm_toy_VA}: length-12 value-added vector
#'   \item \code{bm_toy_X}: length-12 gross output vector
#'   \item \code{bm_toy_countries}: character vector of length 4
#'   \item \code{bm_toy_sectors}: character vector of length 3
#' }
#'
#' The ordering of industries is
#' (China P,M,S; India P,M,S; Japan P,M,S; ROW P,M,S).
#'
#' @format
#' \describe{
#'   \item{bm_toy_Z}{numeric matrix \code{12 x 12}}
#'   \item{bm_toy_Y}{numeric matrix \code{12 x 4}}
#'   \item{bm_toy_VA}{numeric vector of length 12}
#'   \item{bm_toy_X}{numeric vector of length 12}
#'   \item{bm_toy_countries}{character vector of length 4}
#'   \item{bm_toy_sectors}{character vector of length 3}
#' }
#'
#' @docType data
#' @name bm_toy_data
#' @aliases bm_toy_Z bm_toy_Y bm_toy_VA bm_toy_X
#' @aliases bm_toy_countries bm_toy_sectors
#' @keywords datasets
NULL
