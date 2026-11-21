#' Find gene orthologs
#'
#' Find cross-species orthologs for a set of genes.
#'
#' @param genes A character vector of gene symbols or NCBI/Ensembl gene IDs
#'   (case-sensitive).
#' @param source Species of `genes`, such as `"mouse"`. See [species()].
#' @param target Species in which to find orthologs.
#' @param synonyms If `TRUE` (default), check synonyms for input genes.
#' @param best If `TRUE` (default), keep only the best matches.
#' @param min_support Minimum number of sources that support each gene pair.
#' @param orthologs_only If `TRUE` (default), keep only matching pairs.
#' @param details If `TRUE`, return the source gene IDs and the support sources.
#'
#' @return A data frame with one row per distinct input and matching gene pair.
#'   The columns are `input`, `source_symbol`, `target_symbol`,
#'   `target_ncbi_id`, and `target_ensembl_id`. With `details = TRUE`, the
#'   columns also include `source_taxon_id`, `source_ncbi_id`,
#'   `source_ensembl_id`, `target_taxon_id`, `support` (the methods or sources
#'   that support the pair), and `support_n` (the number of these).
#'
#' @details
#' The function uses each distinct input once.
#' Different inputs for the same gene each get their own rows.
#' A valid gene can have no known orthologs in the data.
#'
#' The source AGR dataset includes orthologs from multiple sources, such as
#' PANTHER, SonicParanoid, and OrthoFinder.
#' For pairs supported by three or more methods, a best match in either
#' direction is sufficient for the AGR stringent set. For pairs supported by
#' two methods, a best match in both directions is required. The set also
#' includes all pairs from ZFIN or HGNC.
#'
#' @references
#' Alliance of Genome Resources Consortium. Updates to the Alliance of Genome
#' Resources central infrastructure. Genetics. 2024.
#' doi: \doi{10.1093/genetics/iyae049}.
#'
#' @export
#'
#' @examples
#' orthologs(genes = "TP53", source = "human", target = "mouse")
#' orthologs(genes = "Trp53", source = "mouse", target = "rat")
#' orthologs(genes = c("Ptprc", "Cd34"), source = "mouse", target = "human")
orthologs <- function(
  genes,
  source,
  target,
  synonyms = TRUE,
  best = TRUE,
  min_support = 1,
  orthologs_only = TRUE,
  details = FALSE
) {
  # Check the arguments
  if (!is.atomic(genes) || !is.null(dim(genes))) {
    cli_abort("`genes` must be a vector.")
  }
  # Allow other types, such as integer NCBI IDs, as text
  inputs <- as.character(genes)
  if (is.double(genes) && !is.object(genes)) {
    # Write whole numbers in decimal form so they match the stored NCBI IDs
    # as.character() writes some whole numbers (500000 as "5e+05")
    whole <- is.finite(genes) & genes == trunc(genes)
    inputs[whole] <- sprintf("%.0f", genes[whole])
  }
  genes <- inputs
  if (!rlang::is_bool(synonyms)) {
    cli_abort("`synonyms` must be TRUE or FALSE.")
  }
  if (!rlang::is_bool(best)) {
    cli_abort("`best` must be TRUE or FALSE.")
  }
  if (
    !rlang::is_scalar_integerish(min_support, finite = TRUE) || min_support < 1
  ) {
    cli_abort("`min_support` must be one positive whole number.")
  }
  if (!rlang::is_bool(orthologs_only)) {
    cli_abort("`orthologs_only` must be TRUE or FALSE.")
  }
  if (!rlang::is_bool(details)) {
    cli_abort("`details` must be TRUE or FALSE.")
  }
  if (is.null(source) || is.null(target)) {
    cli_abort("`source` and `target` must each name one species.")
  }
  # Convert species names to taxon IDs
  source_taxon_id <- unique(species(source)$taxon_id)
  target_taxon_id <- unique(species(target)$taxon_id)
  if (source_taxon_id == target_taxon_id) {
    cli_abort("`source` and `target` must name different species.")
  }
  if (anyNA(genes)) {
    cli_warn("Missing values in `genes` were removed.")
    genes <- genes[!is.na(genes)]
  }
  # Search each gene once
  out <- lookup_orthologs(
    unique(genes),
    source_taxon_id,
    target_taxon_id,
    min_support,
    best,
    orthologs_df,
    genes_df,
    synonyms,
    orthologs_only
  )
  if (details) {
    return(out)
  }
  out[c(
    "input",
    "source_symbol",
    "target_symbol",
    "target_ncbi_id",
    "target_ensembl_id"
  )]
}

# Find the orthologs of `genes` from `source` to `target` (both taxon IDs)
# orthologs() checks and cleans the inputs before this step
lookup_orthologs <- function(
  genes,
  source,
  target,
  min_support,
  best,
  pairs,
  gene_table,
  synonyms = FALSE,
  orthologs_only = TRUE
) {
  # Gene 1 is the source only when the source taxon ID is the smaller one
  s <- if (source < target) "g1" else "g2"
  t <- if (source < target) "g2" else "g1"
  pairs <- pairs[
    pairs[[paste0(s, "_taxon_id")]] == source &
      pairs[[paste0(t, "_taxon_id")]] == target,
  ]
  # Rename the stored gene 1 and gene 2 columns to source and target
  pairs <- data.frame(
    source_ncbi_id = pairs[[paste0(s, "_ncbi_id")]],
    target_ncbi_id = pairs[[paste0(t, "_ncbi_id")]],
    support = pairs$algorithms,
    support_n = pairs$algorithms_n,
    best_match = pairs[[paste0(s, "_best_match")]]
  )
  # Match inputs to source genes by NCBI ID, symbol, Ensembl ID, or synonym
  hits <- gene_table[gene_table$ncbi_id %in% pairs$source_ncbi_id, ]
  syn <- synonyms_df[synonyms & synonyms_df$ncbi_id %in% hits$ncbi_id, ]
  matches <- data.frame(
    input = c(
      as.character(hits$ncbi_id),
      hits$symbol,
      hits$ensembl_id,
      syn$synonym
    ),
    source_ncbi_id = c(rep(hits$ncbi_id, 3), syn$ncbi_id)
  )
  matches <- matches[matches$input %in% genes, ]
  # Give each input the pairs of its gene
  out <- dplyr::inner_join(
    dplyr::distinct(matches),
    pairs,
    by = "source_ncbi_id",
    relationship = "many-to-many"
  )
  # Apply the filters
  out <- dplyr::filter(out, .data$support_n >= min_support)
  if (best) {
    out <- dplyr::filter(out, .data$best_match)
  }
  # Inputs without pairs are left out silently, unless no input has pairs
  if (length(genes) && nrow(out) == 0L) {
    cli_warn("No pairs found for any input gene.")
  }
  # Add symbols and Ensembl IDs
  # NCBI IDs are integers in the saved tables and character values in results
  src <- match(out$source_ncbi_id, gene_table$ncbi_id)
  tgt <- match(out$target_ncbi_id, gene_table$ncbi_id)
  out <- data.frame(
    input = out$input,
    source_taxon_id = rep(source, nrow(out)),
    source_symbol = gene_table$symbol[src],
    source_ncbi_id = as.character(out$source_ncbi_id),
    source_ensembl_id = gene_table$ensembl_id[src],
    target_taxon_id = rep(target, nrow(out)),
    target_symbol = gene_table$symbol[tgt],
    target_ncbi_id = as.character(out$target_ncbi_id),
    target_ensembl_id = gene_table$ensembl_id[tgt],
    # The saved table stores the algorithms as a factor to save space
    support = as.character(out$support),
    support_n = out$support_n
  )
  if (!orthologs_only) {
    # Give inputs without pairs one row, with NA in the other columns
    out <- dplyr::bind_rows(
      out,
      data.frame(input = setdiff(genes, out$input))
    )
  }
  # Sort by input order, then by NCBI IDs
  dplyr::arrange(
    out,
    match(.data$input, genes),
    .data$source_ncbi_id,
    .data$target_ncbi_id
  )
}
