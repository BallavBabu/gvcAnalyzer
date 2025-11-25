# R/bm_indices.R

#' Resolve a country (name or index) to numeric index 1..G
#' @export
bm_country_id <- function(io, country) {
  if (!inherits(io, "bm_io")) stop("io must be a 'bm_io' object")

  # numeric
  if (is.numeric(country)) {
    idx <- as.integer(country)
    if (length(idx) != 1L || is.na(idx) || idx < 1L || idx > io$G) {
      stop("Country index out of range 1..", io$G)
    }
    return(idx)
  }

  # character
  if (is.character(country)) {
    if (length(country) != 1L) stop("country must be length-1")
    idx <- io$country_codes[[country]]
    if (is.null(idx)) {
      stop("Unknown country code/name: ", country,
           ". Available: ", paste(io$countries, collapse = ", "))
    }
    return(as.integer(idx))
  }

  stop("country must be numeric index or character code/name")
}

#' Resolve a sector (name or index) to numeric index 1..N
#' @export
bm_sector_id <- function(io, sector) {
  if (!inherits(io, "bm_io")) stop("io must be a 'bm_io' object")

  if (is.numeric(sector)) {
    idx <- as.integer(sector)
    if (length(idx) != 1L || is.na(idx) || idx < 1L || idx > io$N) {
      stop("Sector index out of range 1..", io$N)
    }
    return(idx)
  }

  if (is.character(sector)) {
    if (length(sector) != 1L) stop("sector must be length-1")
    idx <- io$sector_codes[[sector]]
    if (is.null(idx)) {
      stop("Unknown sector code/name: ", sector,
           ". Available: ", paste(io$sectors, collapse = ", "))
    }
    return(as.integer(idx))
  }

  stop("sector must be numeric index or character code/name")
}
