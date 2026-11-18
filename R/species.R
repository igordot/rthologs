#' List available species
#'
#' List all the accepted names for the species in the package data.
#'
#' @param x A species name or NCBI taxon ID.
#'   Use `NULL` (default) to list all species.
#'
#' @return A data frame with `taxon_id`, `name`, and `name_class`.
#'
#' @references
#' Schoch CL, Ciufo S, Domrachev M, et al. NCBI Taxonomy: a comprehensive
#' update on curation, resources and tools. Database (Oxford).
#' 2020. doi: \doi{10.1093/database/baaa062}.
#'
#' @export
#'
#' @examples
#' species()
#' species("Mus musculus")
#' species("rat")
#' species(10090)
species <- function(x = NULL) {
  # Without a species, list all species
  if (is.null(x)) {
    return(species_df)
  }
  # A taxon ID can be a number or a character value
  if (rlang::is_scalar_integerish(x, finite = TRUE)) {
    x <- as.character(as.integer(x))
  }
  if (!rlang::is_string(x)) {
    cli_abort("`x` must be one species name or taxon ID.")
  }
  # Match names without regard to case or extra spaces
  key <- tolower(trimws(x))
  ids <- dplyr::filter(
    species_df,
    tolower(.data$name) == key | .data$taxon_id == key
  )
  ids <- unique(ids$taxon_id)
  if (length(ids) == 0L) {
    cli_abort(
      "Unknown species \"{x}\". Use `species()` to list available names."
    )
  }
  # Return all names of the species
  dplyr::filter(species_df, .data$taxon_id == ids)
}
