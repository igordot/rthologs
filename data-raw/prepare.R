# Prepare orthologs, species, gene identifiers

# Load libraries ----

library(dplyr)
library(tibble)
library(stringr)
library(readr)
library(tidyr)
library(janitor)

# Download source data ----

# Create a directory to store the downloaded files
source_dir <- "downloads"
dir.create(source_dir, showWarnings = FALSE)

# AGR orthologs are the primary source of ortholog pairs
# Using 9.1.0 because the 9.2.0 file includes some missing reverse pairs
agr_orthologs_url <- "https://download.alliancegenome.org/9.1.0/downloads/ORTHOLOGY-ALLIANCE_TSV_COMBINED.tsv.gz"
agr_orthologs_tsv <- file.path(source_dir, "ORTHOLOGY-ALLIANCE.tsv.gz")

# AGR gene records include the primary IDs and cross-references
agr_genes_url <- "https://download.alliancegenome.org/9.1.0/downloads/GENE_TSV_COMBINED.tsv.gz"
agr_genes_tsv <- file.path(source_dir, "GENE.tsv.gz")

# NCBI taxonomy names are used to get species names for the taxon IDs
taxonomy_url <- "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip"
taxonomy_zip <- file.path(source_dir, "taxdmp.zip")

# HGNC is the nomenclature authority for human genes
hgnc_url <- "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
hgnc_tsv <- file.path(source_dir, "hgnc.tsv")

# Allow time for the large downloads
options(timeout = 600)

# Download a file unless it exists (temp name to avoid interrupted downloads)
download_once <- function(url, path) {
  if (!file.exists(path)) {
    part <- str_c(path, ".part")
    download.file(url, part, mode = "wb")
    stopifnot(file.rename(part, path))
  }
  # gzip and zip files store a checksum of their contents
  status <- if (str_ends(path, "\\.zip")) {
    system2("unzip", c("-tqq", shQuote(path)))
  } else if (str_ends(path, "\\.gz")) {
    system2("gzip", c("-t", shQuote(path)))
  } else {
    0
  }
  if (status != 0) {
    stop("Corrupt file, delete it and run again: ", path, call. = FALSE)
  }
}

download_once(agr_orthologs_url, agr_orthologs_tsv)
download_once(agr_genes_url, agr_genes_tsv)
download_once(taxonomy_url, taxonomy_zip)
download_once(hgnc_url, hgnc_tsv)

unzip(taxonomy_zip, files = "names.dmp", exdir = source_dir)

# Species for NCBI and Ensembl imports
species_config <- tribble(
  ~taxon_id , ~ncbi_file                                     , ~ensembl_file                                                         ,
    9606L   , "Mammalia/Homo_sapiens"                        , "homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz"                         ,
   10090L   , "Mammalia/Mus_musculus"                        , "mus_musculus/Mus_musculus.GRCm39.116.gtf.gz"                         ,
   10116L   , "Mammalia/Rattus_norvegicus"                   , "rattus_norvegicus/Rattus_norvegicus.GRCr8.116.gtf.gz"                ,
    7955L   , "Non-mammalian_vertebrates/Danio_rerio"        , "danio_rerio/Danio_rerio.GRCz11.116.gtf.gz"                           ,
    7227L   , "Invertebrates/Drosophila_melanogaster"        , "drosophila_melanogaster/Drosophila_melanogaster.BDGP6.54.63.gtf.gz"  ,
    6239L   , "Invertebrates/Caenorhabditis_elegans"         , "caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.63.gtf.gz"    ,
  559292L   , "Fungi/Saccharomyces_cerevisiae"               , "saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.63.gtf.gz" ,
    8364L   , "Non-mammalian_vertebrates/Xenopus_tropicalis" , "xenopus_tropicalis/Xenopus_tropicalis.UCB_Xtro_10.0.116.gtf.gz"
)
stopifnot(!anyDuplicated(species_config$taxon_id))

# Download NCBI gene records for each species
for (i in seq_len(nrow(species_config))) {
  config <- species_config[i, ]
  download_once(
    str_c(
      "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/",
      config$ncbi_file,
      ".gene_info.gz"
    ),
    file.path(source_dir, str_c("ncbi-", config$taxon_id, ".tsv.gz"))
  )
}

# Download Ensembl GTF files from release 116, the last legacy Ensembl release
# xref.tsv.gz from the new Ensembl FTP site had issues with rat IDs and symbols
for (i in seq_len(nrow(species_config))) {
  config <- species_config[i, ]
  download_once(
    str_c("https://ftp.ensembl.org/pub/release-116/gtf/", config$ensembl_file),
    file.path(source_dir, str_c("ensembl-", config$taxon_id, ".gtf.gz"))
  )
}

# Import source data ----

# Import AGR orthologs
orthologs_all_df <- read_tsv(
  agr_orthologs_tsv,
  skip = 15,
  show_col_types = FALSE
)
stopifnot(nrow(problems(orthologs_all_df)) == 0)
stopifnot(names(orthologs_all_df)[1] == "Gene1ID")
stopifnot(nrow(orthologs_all_df) > 9e5)
glimpse(orthologs_all_df)

# Import AGR gene records
agr_genes_all_df <- read_tsv(
  agr_genes_tsv,
  skip = 14,
  show_col_types = FALSE
)
stopifnot(nrow(problems(agr_genes_all_df)) == 0)
stopifnot(names(agr_genes_all_df)[1] == "Taxon")
stopifnot(nrow(agr_genes_all_df) > 1.6e6)
glimpse(agr_genes_all_df)

# Import NCBI gene records
ncbi_genes_all_df <- tibble()
for (i in seq_len(nrow(species_config))) {
  config <- species_config[i, ]
  ncbi_genes <- withCallingHandlers(
    data.table::fread(
      file = file.path(source_dir, str_c("ncbi-", config$taxon_id, ".tsv.gz")),
      sep = "\t",
      header = TRUE,
      select = c(
        "#tax_id",
        "GeneID",
        "Symbol",
        "dbXrefs",
        "Symbol_from_nomenclature_authority"
      ),
      # Read all columns as character values to prevent type guesses
      colClasses = "character",
      quote = "",
      na.strings = NULL,
      data.table = FALSE,
      nThread = 4
    ),
    warning = function(w) stop(conditionMessage(w), call. = FALSE)
  ) |>
    as_tibble()
  stopifnot(config$taxon_id %in% ncbi_genes$`#tax_id`)
  ncbi_genes_all_df <- bind_rows(ncbi_genes_all_df, ncbi_genes)
}
stopifnot(nrow(ncbi_genes_all_df) > 5e5)
glimpse(ncbi_genes_all_df)
count(ncbi_genes_all_df, `#tax_id`)

# Import Ensembl gene rows from the GTF files
ensembl_genes_all_df <- tibble()
for (i in seq_len(nrow(species_config))) {
  config <- species_config[i, ]
  gtf_file <- file.path(
    source_dir,
    str_c("ensembl-", config$taxon_id, ".gtf.gz")
  )
  gtf <- withCallingHandlers(
    data.table::fread(
      # The system gzip is faster than the decompression in fread()
      cmd = str_c("gzip -dc ", shQuote(gtf_file)),
      sep = "\t",
      header = FALSE,
      # Start at the first feature row, after the "#!" header lines
      skip = "gene_id",
      select = c(3, 9),
      col.names = c("feature", "attributes"),
      quote = "",
      na.strings = NULL,
      data.table = FALSE,
      nThread = 4
    ),
    warning = function(w) stop(conditionMessage(w), call. = FALSE)
  )
  # Attributes of the gene rows (not transcripts or exons)
  attributes <- gtf$attributes[gtf$feature == "gene"]
  ensembl_genes <- tibble(
    taxon_id = config$taxon_id,
    # Try to match whole keys, not keys that end in gene_id or gene_name
    gene_id = str_match(attributes, '(?:^|; )gene_id "([^"]+)"')[, 2],
    gene_name = str_match(attributes, '(?:^|; )gene_name "([^"]+)"')[, 2]
  )
  stopifnot(nrow(ensembl_genes) > 5000)
  stopifnot(!anyNA(ensembl_genes$gene_id))
  stopifnot(!anyDuplicated(ensembl_genes$gene_id))
  ensembl_genes_all_df <- bind_rows(ensembl_genes_all_df, ensembl_genes)
}
glimpse(ensembl_genes_all_df)
count(ensembl_genes_all_df, taxon_id)

# Import HGNC gene records
hgnc_genes_all_df <- read_tsv(
  hgnc_tsv,
  col_types = cols(.default = col_character())
)
stopifnot(nrow(problems(hgnc_genes_all_df)) == 0)
stopifnot(c("hgnc_id", "ensembl_gene_id") %in% names(hgnc_genes_all_df))
stopifnot(nrow(hgnc_genes_all_df) > 4e4)
stopifnot(!anyDuplicated(hgnc_genes_all_df$hgnc_id))
glimpse(hgnc_genes_all_df)

# Import NCBI taxonomy names
taxonomy_df <- read_delim(
  file.path(source_dir, "names.dmp"),
  delim = "|",
  trim_ws = TRUE,
  col_names = FALSE,
  show_col_types = FALSE
)
stopifnot(nrow(problems(taxonomy_df)) == 0)
stopifnot(nrow(taxonomy_df) > 5e6)
glimpse(taxonomy_df)

# Prepare orthologs table ----

# Simplify column names
orthologs_df <- clean_names(orthologs_all_df)
names(orthologs_df) <- str_replace(names(orthologs_df), "^gene([12])", "g\\1_")
names(orthologs_df) <- str_replace(
  names(orthologs_df),
  "^g([12])_id$",
  "g\\1_agr_id"
)
names(orthologs_df) <- str_replace(
  names(orthologs_df),
  "_species_taxon_id$",
  "_taxon_id"
)

# Convert taxon IDs to integers
orthologs_df <- mutate(
  orthologs_df,
  g1_taxon_id = str_remove(g1_taxon_id, "NCBITaxon:"),
  g1_taxon_id = as.integer(g1_taxon_id),
  g2_taxon_id = str_remove(g2_taxon_id, "NCBITaxon:"),
  g2_taxon_id = as.integer(g2_taxon_id),
  algorithms_match = as.integer(algorithms_match)
)
stopifnot(all(orthologs_df$g1_taxon_id %in% orthologs_df$g2_taxon_id))
stopifnot(all(orthologs_df$g2_taxon_id %in% orthologs_df$g1_taxon_id))
stopifnot(!anyNA(orthologs_df[c("g1_taxon_id", "g2_taxon_id")]))

# Select the necessary columns
relevant_columns <- c(
  "g1_agr_id",
  "g1_symbol",
  "g1_taxon_id",
  "g1_species_name",
  "g2_agr_id",
  "g2_symbol",
  "g2_taxon_id",
  "g2_species_name",
  "algorithms",
  "algorithms_match",
  "out_of_algorithms",
  "is_best_score"
)
stopifnot(all(relevant_columns %in% names(orthologs_df)))
orthologs_df <- orthologs_df[relevant_columns]
stopifnot(!anyNA(orthologs_df))
stopifnot(all(nzchar(as.matrix(orthologs_df))))
glimpse(orthologs_df)

# Remove species with few partners (some species have limited pairs)
tax_pairs <- unique(data.frame(
  first = pmin(orthologs_df$g1_taxon_id, orthologs_df$g2_taxon_id),
  second = pmax(orthologs_df$g1_taxon_id, orthologs_df$g2_taxon_id)
))
partners <- table(c(tax_pairs$first, tax_pairs$second))
excluded <- names(partners)[partners < 5]
if (length(excluded)) {
  message("Exclude taxa: ", toString(excluded))
}
tax_pairs <- tax_pairs[
  !tax_pairs$first %in% excluded & !tax_pairs$second %in% excluded,
]
partners <- table(c(tax_pairs$first, tax_pairs$second))
stopifnot(length(unique(partners)) == 1)

# Remaining species must all have the same number of partners
orthologs_df <- filter(
  orthologs_df,
  !g1_taxon_id %in% excluded,
  !g2_taxon_id %in% excluded
)

# Convert the best score column to logical values
stopifnot(all(orthologs_df$is_best_score %in% c("Yes", "No", "Yes_Adjusted")))
orthologs_df$is_best_score <- orthologs_df$is_best_score %in%
  c("Yes", "Yes_Adjusted")
glimpse(orthologs_df)

# Keep species names and directional flags before merging reverse pairs
ortholog_species <- distinct(
  orthologs_df,
  taxon_id = g1_taxon_id,
  name = g1_species_name,
  name_class = "species name"
)
best_score <- orthologs_df$is_best_score

# Store each pair with the smaller taxon ID first
stopifnot(all(orthologs_df$g1_taxon_id != orthologs_df$g2_taxon_id))
swap <- orthologs_df$g1_taxon_id > orthologs_df$g2_taxon_id
orthologs_df <- data.frame(
  g1_taxon_id = pmin(
    orthologs_df$g1_taxon_id,
    orthologs_df$g2_taxon_id
  ),
  g1_agr_id = if_else(swap, orthologs_df$g2_agr_id, orthologs_df$g1_agr_id),
  g1_symbol = if_else(swap, orthologs_df$g2_symbol, orthologs_df$g1_symbol),
  g2_taxon_id = pmax(
    orthologs_df$g1_taxon_id,
    orthologs_df$g2_taxon_id
  ),
  g2_agr_id = if_else(swap, orthologs_df$g1_agr_id, orthologs_df$g2_agr_id),
  g2_symbol = if_else(swap, orthologs_df$g1_symbol, orthologs_df$g2_symbol),
  algorithms = orthologs_df$algorithms,
  algorithms_n = orthologs_df$algorithms_match,
  algorithms_prop = round(
    orthologs_df$algorithms_match / orthologs_df$out_of_algorithms,
    3
  )
)
glimpse(orthologs_df)

# Check that each gene pair appears twice (once in each direction)
keys <- c("g1_taxon_id", "g1_agr_id", "g2_taxon_id", "g2_agr_id")
pair_key <- do.call(paste, orthologs_df[keys])
stopifnot(all(table(pair_key) == 2))

# Check that algorithm totals agree in both directions
reverse_index <- match(pair_key[!swap], pair_key[swap])
stopifnot(!anyNA(reverse_index))
stopifnot(identical(
  orthologs_df$algorithms[!swap],
  orthologs_df$algorithms[swap][reverse_index]
))

# Remove duplicate gene pairs
orthologs_df <- distinct(orthologs_df)
stopifnot(!anyDuplicated(orthologs_df[keys]))
stopifnot(length(swap) == nrow(orthologs_df) * 2)

# Most gene pairs should be retained
stopifnot(
  mean(
    unique(paste(
      pmin(orthologs_all_df$Gene1ID, orthologs_all_df$Gene2ID),
      pmax(orthologs_all_df$Gene1ID, orthologs_all_df$Gene2ID)
    )) %in%
      paste(
        pmin(orthologs_df$g1_agr_id, orthologs_df$g2_agr_id),
        pmax(orthologs_df$g1_agr_id, orthologs_df$g2_agr_id)
      )
  ) >=
    0.9
)

# Take each best-match flag from its own direction (they can disagree)
forward_index <- match(do.call(paste, orthologs_df[keys]), pair_key[!swap])
backward_index <- match(do.call(paste, orthologs_df[keys]), pair_key[swap])
stopifnot(!anyNA(forward_index))
stopifnot(!anyNA(backward_index))
orthologs_df$g1_best_match <- best_score[!swap][forward_index]
orthologs_df$g2_best_match <- best_score[swap][backward_index]
rownames(orthologs_df) <- NULL
glimpse(orthologs_df)

# Prepare species table ----

# Get the relevant taxon IDs from the orthologs table
tax_ids <- unique(c(orthologs_df$g1_taxon_id, orthologs_df$g2_taxon_id))
stopifnot(setequal(tax_ids, species_config$taxon_id))

# Keep only the relevant taxon IDs in the taxonomy table
stopifnot(all(tax_ids %in% taxonomy_df$X1))
taxonomy_df <- filter(taxonomy_df, X1 %in% tax_ids)

# Create a species data frame with taxon IDs and all name variants
species_df <- select(
  taxonomy_df,
  taxon_id = X1,
  name = X2,
  name_class = X4
)
species_df$taxon_id <- as.integer(species_df$taxon_id)
species_df <- filter(species_df, taxon_id %in% tax_ids)
species_df <- filter(species_df, str_detect(name_class, "name"))
species_df <- bind_rows(species_df, ortholog_species)
species_df <- distinct(species_df)
stopifnot(all(species_df$taxon_id %in% tax_ids))
glimpse(species_df)

# Prepare AGR genes table ----

# Simplify AGR genes table
agr_genes_df <- agr_genes_all_df |>
  clean_names() |>
  rename(
    taxon_id = taxon,
    agr_id = gene_id,
    symbol = gene_symbol,
    synonyms = gene_synonyms,
    cross_references = gene_cross_references
  ) |>
  mutate(taxon_id = as.integer(str_remove(taxon_id, "^NCBITaxon:")))
stopifnot(all(tax_ids %in% agr_genes_df$taxon_id))
stopifnot(sum(agr_genes_df$taxon_id %in% tax_ids) > 1.6e6)

# AGR IDs used by the ortholog pairs
agr_gene_ids <- sort(unique(c(orthologs_df$g1_agr_id, orthologs_df$g2_agr_id)))
stopifnot(
  !anyDuplicated(agr_genes_df$agr_id[agr_genes_df$taxon_id %in% tax_ids])
)
stopifnot(all(agr_gene_ids %in% agr_genes_df$agr_id))

# Keep only the ortholog genes and the necessary columns
agr_genes_df <- agr_genes_df |>
  filter(agr_id %in% agr_gene_ids) |>
  select(taxon_id, agr_id, symbol, synonyms, cross_references)
stopifnot(!anyDuplicated(agr_genes_df$agr_id))

# Expand cross-references to one identifier per row
agr_xrefs <- agr_genes_df |>
  transmute(agr_id, cross_reference = cross_references) |>
  separate_longer_delim(cross_reference, delim = "|") |>
  filter(!is.na(cross_reference), str_detect(cross_reference, ":")) |>
  separate_wider_delim(
    cross_reference,
    delim = ":",
    names = c("database", "external_id"),
    too_many = "merge"
  ) |>
  filter(!is.na(external_id), !external_id %in% c("", "-"))
glimpse(agr_xrefs)
count(agr_xrefs, database)

# Prepare NCBI genes table ----

# Simplify NCBI genes table and match the column names to other tables
ncbi_genes_df <- ncbi_genes_all_df |>
  clean_names() |>
  rename(taxon_id = number_tax_id) |>
  rename(
    ncbi_id = gene_id,
    authority_symbol = symbol_from_nomenclature_authority,
    cross_references = db_xrefs
  ) |>
  mutate(taxon_id = as.integer(taxon_id), ncbi_id = as.integer(ncbi_id))
stopifnot(!anyNA(ncbi_genes_df$ncbi_id))
stopifnot(all(tax_ids %in% ncbi_genes_df$taxon_id))

# Keep only the relevant taxa
ncbi_genes_df <- filter(ncbi_genes_df, taxon_id %in% tax_ids)
stopifnot(nrow(ncbi_genes_df) > 5e5)

# Expand cross-references to one identifier per row
ncbi_all_xrefs <- ncbi_genes_df |>
  transmute(taxon_id, ncbi_id, cross_reference = cross_references) |>
  separate_longer_delim(cross_reference, delim = "|") |>
  filter(!is.na(cross_reference), str_detect(cross_reference, ":")) |>
  separate_wider_delim(
    cross_reference,
    delim = ":",
    names = c("database", "external_id"),
    too_many = "merge"
  )

# Cross-references with database names normalized to match AGR
ncbi_all_xrefs <- mutate(
  ncbi_all_xrefs,
  database = replace_values(database, "FLYBASE" ~ "FB", "WormBase" ~ "WB")
)
glimpse(ncbi_all_xrefs)
count(ncbi_all_xrefs, database)

# Match NCBI cross-references to AGR primary IDs
ncbi_agr_xrefs <- ncbi_all_xrefs |>
  mutate(
    external_id = case_when(
      database == "HGNC" ~ str_remove(external_id, "^HGNC:"),
      database == "MGI" ~ str_remove(external_id, "^MGI:"),
      .default = external_id
    ),
    agr_id = paste(database, external_id, sep = ":")
  ) |>
  filter(agr_id %in% agr_gene_ids)
stopifnot(all(tax_ids %in% ncbi_agr_xrefs$taxon_id))
glimpse(ncbi_agr_xrefs)
count(ncbi_agr_xrefs, database)

# Select NCBI cross-references to Ensembl IDs (not all species)
# Ensembl uses the species database gene ID for FlyBase, WormBase, SGD
ncbi_ensembl_xrefs <- ncbi_all_xrefs |>
  filter(database %in% c("Ensembl", "EnsemblRapid")) |>
  select(taxon_id, ncbi_id, ensembl_id = external_id)
glimpse(ncbi_ensembl_xrefs)
count(ncbi_ensembl_xrefs, taxon_id)

# Prepare Ensembl genes table ----

# Use the gene ID when there is no gene name
ensembl_genes_df <- transmute(
  ensembl_genes_all_df,
  taxon_id,
  ensembl_id = gene_id,
  symbol = coalesce(gene_name, gene_id)
)
glimpse(ensembl_genes_df)
count(ensembl_genes_df, taxon_id)

# Prepare ortholog genes ----

# NCBI IDs and Ensembl IDs use two separate sets of steps
# Each set starts with the single IDs from AGR and keeps them unchanged
# Later steps use only the genes that remain without a selected ID

# Get the AGR cross-references to one database as a list for each ortholog gene
agr_id_lists <- function(database_name) {
  xrefs <- filter(agr_xrefs, database == database_name)
  ids <- split(xrefs$external_id, xrefs$agr_id)[agr_genes_df$agr_id]
  # Genes without a cross-reference get an empty vector
  unname(lapply(ids, function(x) unique(as.character(x))))
}

# Add NCBI IDs step 1: single NCBI IDs from AGR ----

# Genes with their symbols and the NCBI IDs listed by AGR
ncbi_pending <- agr_genes_df |>
  select(taxon_id, agr_id, symbol) |>
  mutate(
    agr_ncbi_ids = lapply(agr_id_lists("NCBI_Gene"), as.integer),
    ncbi_id = NA_integer_,
    ncbi_source = NA_character_
  )
stopifnot(!anyNA(unlist(ncbi_pending$agr_ncbi_ids)))

# Accept genes with exactly one AGR NCBI ID
single <- lengths(ncbi_pending$agr_ncbi_ids) == 1
ncbi_pending$ncbi_id[single] <- unlist(ncbi_pending$agr_ncbi_ids[single])
ncbi_pending$ncbi_source[single] <- "AGR"

# Genes with an NCBI ID selected in step 1
ncbi_step1 <- filter(ncbi_pending, !is.na(ncbi_id))
nrow(ncbi_step1)

# Genes still without an NCBI ID (AGR lists zero or multiple IDs)
ncbi_pending <- filter(ncbi_pending, is.na(ncbi_id))
count(ncbi_pending, taxon_id, agr_ncbi_ids = lengths(agr_ncbi_ids))
stopifnot(nrow(ncbi_pending) < 500)

# Add NCBI IDs step 2: choose among multiple AGR NCBI IDs ----

# Keep the AGR NCBI ID whose NCBI symbol matches the AGR symbol
ncbi_matches <- ncbi_pending |>
  select(agr_id, taxon_id, symbol, ncbi_id = agr_ncbi_ids) |>
  unnest_longer(ncbi_id, ptype = integer()) |>
  left_join(
    select(
      ncbi_genes_df,
      taxon_id,
      ncbi_id,
      ncbi_symbol = symbol,
      authority_symbol
    ),
    by = c("taxon_id", "ncbi_id")
  ) |>
  filter(symbol == ncbi_symbol | symbol == authority_symbol) |>
  transmute(agr_id, ncbi_id, ncbi_source = "AGR")
stopifnot(!anyNA(ncbi_matches))
stopifnot(!anyDuplicated(ncbi_matches$agr_id))
ncbi_pending <- rows_update(ncbi_pending, ncbi_matches, by = "agr_id")

# Genes with an NCBI ID selected in step 2
ncbi_step2 <- filter(ncbi_pending, !is.na(ncbi_id))
nrow(ncbi_step2)

# Genes still without an NCBI ID
ncbi_pending <- filter(ncbi_pending, is.na(ncbi_id))
nrow(ncbi_pending)

# Add NCBI IDs step 3: add NCBI IDs of NCBI genes that list the AGR ID ----

# Use NCBI genes that link back to AGR genes with no AGR NCBI ID
ncbi_matches <- ncbi_pending |>
  filter(lengths(agr_ncbi_ids) == 0) |>
  select(agr_id, taxon_id) |>
  inner_join(
    distinct(ncbi_agr_xrefs, taxon_id, agr_id, ncbi_id),
    by = c("taxon_id", "agr_id")
  ) |>
  transmute(agr_id, ncbi_id, ncbi_source = "NCBI")
stopifnot(!anyNA(ncbi_matches))
stopifnot(!anyDuplicated(ncbi_matches$agr_id))
ncbi_pending <- rows_update(ncbi_pending, ncbi_matches, by = "agr_id")

# Genes with an NCBI ID selected in step 3
ncbi_step3 <- filter(ncbi_pending, !is.na(ncbi_id))
nrow(ncbi_step3)

# Genes still without an NCBI ID
ncbi_pending <- filter(ncbi_pending, is.na(ncbi_id))
ncbi_pending

# Few genes should remain without an NCBI ID
stopifnot(nrow(ncbi_pending) < 10)

# Combine the NCBI steps ----

# All ortholog genes with the selected NCBI ID and the step that selected it
genes_ncbi <- bind_rows(
  "1 one in AGR" = ncbi_step1,
  "2 multiple in AGR" = ncbi_step2,
  "3 none in AGR" = ncbi_step3,
  "unresolved" = ncbi_pending,
  .id = "ncbi_step"
)
stopifnot(setequal(agr_gene_ids, genes_ncbi$agr_id))
stopifnot(!anyDuplicated(genes_ncbi$agr_id))
stopifnot(identical(is.na(genes_ncbi$ncbi_id), is.na(genes_ncbi$ncbi_source)))

# Count NCBI ID selections by step and species
ncbi_steps_summary <- genes_ncbi |>
  count(taxon_id, ncbi_step) |>
  pivot_wider(names_from = taxon_id, values_from = n, values_fill = 0) |>
  arrange(ncbi_step)
ncbi_steps_summary

# Add Ensembl IDs step 1: single Ensembl IDs from AGR ----

# Genes with their symbols, selected NCBI IDs, and Ensembl IDs listed by AGR
ensembl_pending <- agr_genes_df |>
  select(taxon_id, agr_id, symbol) |>
  mutate(
    agr_ensembl_ids = agr_id_lists("ENSEMBL"),
    ensembl_candidates = agr_ensembl_ids,
    ensembl_id = NA_character_,
    ensembl_source = NA_character_
  ) |>
  left_join(
    select(genes_ncbi, agr_id, ncbi_id),
    by = "agr_id",
    relationship = "one-to-one"
  )

# Select the candidate of pending genes with exactly one Ensembl candidate
accept_single <- function(genes, source) {
  single <- is.na(genes$ensembl_id) & lengths(genes$ensembl_candidates) == 1
  genes$ensembl_id[single] <- unlist(genes$ensembl_candidates[single])
  genes$ensembl_source[single] <- source
  genes
}

# Accept genes with exactly one AGR Ensembl ID
ensembl_pending <- accept_single(ensembl_pending, "AGR")

# Genes with an Ensembl ID selected in step 1
ensembl_step1 <- filter(ensembl_pending, !is.na(ensembl_id))
nrow(ensembl_step1)

# Genes still without an Ensembl ID (AGR lists zero or multiple IDs)
ensembl_pending <- filter(ensembl_pending, is.na(ensembl_id))
count(ensembl_pending, taxon_id)

# Add Ensembl IDs step 2: add the Ensembl IDs from HGNC lists ----

# HGNC lists at most one Ensembl ID for each human gene
# AGR human gene IDs are HGNC IDs
hgnc_ensembl <- hgnc_genes_all_df |>
  filter(!is.na(ensembl_gene_id)) |>
  select(agr_id = hgnc_id, ensembl_id = ensembl_gene_id)
stopifnot(all(str_detect(hgnc_ensembl$ensembl_id, "^ENSG[0-9]{11}$")))
stopifnot(all(str_detect(
  ensembl_pending$agr_id[ensembl_pending$taxon_id == 9606],
  "^HGNC:[0-9]+$"
)))

# Use the HGNC Ensembl ID for pending human genes
hgnc_matches <- hgnc_ensembl |>
  filter(agr_id %in% ensembl_pending$agr_id) |>
  mutate(
    ensembl_source = "HGNC",
    ensembl_candidates = as.list(ensembl_id)
  )
ensembl_pending <- rows_update(ensembl_pending, hgnc_matches, by = "agr_id")

# Genes with an Ensembl ID selected in step 2
ensembl_step2 <- filter(ensembl_pending, !is.na(ensembl_id))
nrow(ensembl_step2)

# Genes still without an Ensembl ID
ensembl_pending <- filter(ensembl_pending, is.na(ensembl_id))
count(ensembl_pending, taxon_id)

# Add Ensembl IDs step 3: keep the AGR Ensembl IDs that Ensembl has ----

# Some AGR Ensembl IDs are not in the Ensembl GTF file of the species
candidate_ids <- unlist(ensembl_pending$ensembl_candidates)
candidate_rows <- rep(
  seq_len(nrow(ensembl_pending)),
  lengths(ensembl_pending$ensembl_candidates)
)
in_ensembl <- paste(
  ensembl_pending$taxon_id[candidate_rows],
  candidate_ids
) %in%
  paste(ensembl_genes_df$taxon_id, ensembl_genes_df$ensembl_id)
ensembl_pending$ensembl_candidates <- unname(split(
  candidate_ids[in_ensembl],
  factor(candidate_rows[in_ensembl], levels = seq_len(nrow(ensembl_pending)))
))

# Accept genes with exactly one AGR Ensembl ID left
ensembl_pending <- accept_single(ensembl_pending, "Ensembl")

# Genes with an Ensembl ID selected in step 3
ensembl_step3 <- filter(ensembl_pending, !is.na(ensembl_id))
nrow(ensembl_step3)
count(ensembl_step3, taxon_id)

# Genes still without an Ensembl ID
ensembl_pending <- filter(ensembl_pending, is.na(ensembl_id))
count(ensembl_pending, taxon_id)

# Add Ensembl IDs step 4: choose among multiple AGR Ensembl IDs by symbol ----

# Sources can write one symbol in different cases (Ypel1 versus YPEL1)
# Find those symbols so that the lookup uses only exact matches for them
agr_symbols <- distinct(agr_genes_df, taxon_id, symbol)
agr_symbol_keys <- str_to_lower(paste(agr_symbols$taxon_id, agr_symbols$symbol))
case_conflicts <- unique(agr_symbol_keys[duplicated(agr_symbol_keys)])
case_conflicts

# The only known conflict is the Xenbase placeholder (XXX versus xxx)
stopifnot(identical(case_conflicts, "8364 xxx"))

# Ensembl IDs for each species and exact symbol
symbol_ids <- split(
  ensembl_genes_df$ensembl_id,
  paste(ensembl_genes_df$taxon_id, ensembl_genes_df$symbol)
)

# Ensembl IDs for each species and symbol when the lookup ignores case
symbol_ids_any_case <- split(
  ensembl_genes_df$ensembl_id,
  str_to_lower(paste(ensembl_genes_df$taxon_id, ensembl_genes_df$symbol))
)

# Find the Ensembl IDs whose Ensembl gene name matches a symbol
find_by_symbol <- function(taxon, symbol) {
  ids <- symbol_ids[[paste(taxon, symbol)]]
  # Ignore case only when there is no exact match and no case conflict
  key <- str_to_lower(paste(taxon, symbol))
  if (is.null(ids) && !key %in% case_conflicts) {
    ids <- symbol_ids_any_case[[key]]
  }
  as.character(ids)
}

# Keep only the candidates that are also in `ids`
# If no candidate is in `ids`, keep all the candidates
restrict <- function(candidates, ids) {
  matched <- intersect(candidates, ids)
  if (length(matched)) matched else candidates
}

# Keep the candidates whose Ensembl gene name matches the AGR symbol
rows <- which(lengths(ensembl_pending$ensembl_candidates) > 1)
ensembl_pending$ensembl_candidates[rows] <- with(
  ensembl_pending[rows, ],
  Map(
    function(taxon, symbol, candidates) {
      restrict(candidates, find_by_symbol(taxon, symbol))
    },
    taxon_id,
    symbol,
    ensembl_candidates
  )
)
ensembl_pending <- accept_single(ensembl_pending, "Ensembl")

# Genes with an Ensembl ID selected in step 4
ensembl_step4 <- filter(ensembl_pending, !is.na(ensembl_id))
nrow(ensembl_step4)
count(ensembl_step4, taxon_id)

# Genes still without an Ensembl ID
ensembl_pending <- filter(ensembl_pending, is.na(ensembl_id))
count(ensembl_pending, taxon_id)

# Add Ensembl IDs step 5: choose among multiple AGR Ensembl IDs by NCBI ----

# Ensembl IDs that NCBI lists for each species and NCBI ID
ncbi_ensembl_ids <- split(
  ncbi_ensembl_xrefs$ensembl_id,
  paste(ncbi_ensembl_xrefs$taxon_id, ncbi_ensembl_xrefs$ncbi_id)
)

# Keep the candidates that NCBI lists for the selected NCBI ID
rows <- which(lengths(ensembl_pending$ensembl_candidates) > 1)
ensembl_pending$ensembl_candidates[rows] <- with(
  ensembl_pending[rows, ],
  Map(
    function(taxon, ncbi_id, candidates) {
      restrict(candidates, ncbi_ensembl_ids[[paste(taxon, ncbi_id)]])
    },
    taxon_id,
    ncbi_id,
    ensembl_candidates
  )
)
ensembl_pending <- accept_single(ensembl_pending, "NCBI")

# Genes with an Ensembl ID selected in step 5
ensembl_step5 <- filter(ensembl_pending, !is.na(ensembl_id))
nrow(ensembl_step5)
count(ensembl_step5, taxon_id)

# Genes still without an Ensembl ID
ensembl_pending <- filter(ensembl_pending, is.na(ensembl_id))
count(ensembl_pending, taxon_id)

# Add Ensembl IDs step 6: add Ensembl IDs that NCBI lists for the NCBI ID ----

# Skip genes for which NCBI lists more than one usable Ensembl ID
ncbi_ensembl_matches <- ensembl_pending |>
  filter(lengths(ensembl_candidates) == 0, !is.na(ncbi_id)) |>
  select(agr_id, taxon_id, ncbi_id) |>
  inner_join(
    ncbi_ensembl_xrefs,
    by = c("taxon_id", "ncbi_id"),
    relationship = "many-to-many"
  ) |>
  filter(!ensembl_id %in% c("", "-")) |>
  distinct(agr_id, ensembl_id) |>
  filter(n() == 1, .by = agr_id) |>
  mutate(
    ensembl_source = "NCBI",
    ensembl_candidates = as.list(ensembl_id)
  )
ensembl_pending <- rows_update(
  ensembl_pending,
  ncbi_ensembl_matches,
  by = "agr_id"
)

# Genes with an Ensembl ID selected in step 6
ensembl_step6 <- filter(ensembl_pending, !is.na(ensembl_id))
nrow(ensembl_step6)
count(ensembl_step6, taxon_id)

# Genes still without an Ensembl ID
ensembl_pending <- filter(ensembl_pending, is.na(ensembl_id))
count(ensembl_pending, taxon_id)

# Add Ensembl IDs step 7: add Ensembl IDs from Ensembl by symbol ----

# AGR uses fly FlyBase IDs, AGR and NCBI do not list them as Ensembl IDs
# Symbols are inconsistent between databases (14-3-3ε versus 14-3-3epsilon)
fly_genes <- filter(ensembl_pending, taxon_id == 7227)
fly_gtf_ids <- ensembl_genes_df$ensembl_id[ensembl_genes_df$taxon_id == 7227]

# AGR and NCBI give no fly Ensembl IDs
stopifnot(nrow(fly_genes) > 0)
stopifnot(all(lengths(fly_genes$agr_ensembl_ids) == 0))
stopifnot(!7227 %in% ncbi_ensembl_xrefs$taxon_id)

# AGR fly IDs are FlyBase gene IDs
stopifnot(all(str_detect(fly_genes$agr_id, "^FB:FBgn[0-9]{7}$")))
fly_flybase_ids <- str_remove(fly_genes$agr_id, "^FB:")
stopifnot(mean(fly_flybase_ids %in% fly_gtf_ids) > 0.99)

# Fly GTF records whose symbol identifies one gene
fly_symbol_ids <- ensembl_genes_df |>
  filter(taxon_id == 7227) |>
  filter(n() == 1, .by = symbol)

# Ensembl ID that the symbol lookup finds for each pending fly gene
fly_symbol_matches <- fly_symbol_ids$ensembl_id[
  match(fly_genes$symbol, fly_symbol_ids$symbol)
]

# Check the symbol lookup
stopifnot(anyNA(fly_symbol_matches))
stopifnot(all(fly_symbol_matches == fly_flybase_ids, na.rm = TRUE))

# Ensembl IDs that earlier steps gave to other genes
taken <- c(
  ensembl_step1$ensembl_id,
  ensembl_step2$ensembl_id,
  ensembl_step3$ensembl_id,
  ensembl_step4$ensembl_id,
  ensembl_step5$ensembl_id,
  ensembl_step6$ensembl_id
)

# NCBI IDs that list each Ensembl ID
linked_ncbi <- split(ncbi_ensembl_xrefs$ncbi_id, ncbi_ensembl_xrefs$ensembl_id)

# Find the Ensembl ID candidates for a gene without a usable AGR Ensembl ID
find_ensembl <- function(taxon, agr_id, symbol, ncbi_id) {
  # Fly exception: use the FlyBase ID (see the checks above)
  if (taxon == 7227) {
    return(intersect(str_remove(agr_id, "^FB:"), fly_gtf_ids))
  }
  ids <- find_by_symbol(taxon, symbol)
  # Skip IDs given to other genes and IDs that NCBI links to other genes
  ncbi_ok <- vapply(
    ids,
    function(id) all(linked_ncbi[[id]] %in% ncbi_id),
    logical(1)
  )
  ids[!ids %in% taken & ncbi_ok]
}

# Find the candidates for each pending gene without candidates
rows <- which(lengths(ensembl_pending$ensembl_candidates) == 0)
ensembl_pending$ensembl_candidates[rows] <- with(
  ensembl_pending[rows, ],
  Map(find_ensembl, taxon_id, agr_id, symbol, ncbi_id)
)

# Accept only the genes of this step
# Genes with multiple AGR candidates from step 5 stay unresolved
single <- rows[lengths(ensembl_pending$ensembl_candidates[rows]) == 1]
ensembl_pending$ensembl_id[single] <- unlist(
  ensembl_pending$ensembl_candidates[single]
)
ensembl_pending$ensembl_source[single] <- "Ensembl"

# Genes with an Ensembl ID selected in step 7
ensembl_step7 <- filter(ensembl_pending, !is.na(ensembl_id))
nrow(ensembl_step7)
count(ensembl_step7, taxon_id)

# Genes still without an Ensembl ID
ensembl_pending <- filter(ensembl_pending, is.na(ensembl_id))
count(ensembl_pending, taxon_id)

# Combine the Ensembl steps ----

# All ortholog genes with the selected Ensembl ID and the step that selected it
genes_ensembl <- bind_rows(
  "1 one in AGR" = ensembl_step1,
  "2 none or multiple in AGR, HGNC" = ensembl_step2,
  "3 multiple in AGR, one in Ensembl" = ensembl_step3,
  "4 multiple in AGR, Ensembl symbol" = ensembl_step4,
  "5 multiple in AGR, NCBI" = ensembl_step5,
  "6 none in AGR, NCBI" = ensembl_step6,
  "7 none in AGR, Ensembl symbol" = ensembl_step7,
  "unresolved" = ensembl_pending,
  .id = "ensembl_step"
)
stopifnot(setequal(agr_gene_ids, genes_ensembl$agr_id))
stopifnot(!anyDuplicated(genes_ensembl$agr_id))
stopifnot(identical(
  is.na(genes_ensembl$ensembl_id),
  is.na(genes_ensembl$ensembl_source)
))

# Count Ensembl ID selections by step and species
ensembl_steps_summary <- genes_ensembl |>
  count(taxon_id, ensembl_step) |>
  pivot_wider(names_from = taxon_id, values_from = n, values_fill = 0) |>
  arrange(ensembl_step)
ensembl_steps_summary

# Combine the NCBI and Ensembl gene tables ----

# Gene records with selected identifiers and their sources
genes_df <- agr_genes_df |>
  select(taxon_id, agr_id, symbol, synonyms) |>
  left_join(
    select(genes_ncbi, agr_id, ncbi_id, ncbi_source),
    by = "agr_id",
    relationship = "one-to-one"
  ) |>
  left_join(
    select(genes_ensembl, agr_id, ensembl_id, ensembl_source),
    by = "agr_id",
    relationship = "one-to-one"
  )
stopifnot(identical(genes_df$agr_id, agr_genes_df$agr_id))

# Prepare gene synonyms ----

# Expand AGR synonyms to one row per gene and synonym
synonyms_df <- tibble(
  taxon_id = genes_df$taxon_id,
  agr_id = genes_df$agr_id,
  synonym = strsplit(genes_df$synonyms, "|", fixed = TRUE)
) |>
  unnest_longer(synonym, ptype = character()) |>
  distinct()
nrow(synonyms_df)

# Remove missing values or placeholders
# Remove whitespace, list separators (, * ;), or ? and &
# (? and & are mostly damaged Greek letters, such as Pk?1 and &agr-PS3)
# Remove digits only because they look like NCBI Gene IDs
# Remove tags, such as p19<sup>ARF</sup>
# Remove worm locus tags, such as CELE_Y74C9A.3
# Remove synonyms longer than 25 characters (the longest main symbol has 18)
synonyms_df <- filter(
  synonyms_df,
  !is.na(synonym),
  nzchar(synonym),
  synonym != "-",
  !str_detect(synonym, "[\\s,*;?\\&]"),
  !str_detect(synonym, "^[0-9]+$"),
  !str_detect(synonym, "<[^>]+>"),
  !str_starts(synonym, fixed("CELE_")),
  str_length(synonym) <= 25
)
nrow(synonyms_df)

# Remove synonyms that match any primary symbol in the same species
synonyms_df <- anti_join(
  synonyms_df,
  transmute(genes_df, taxon_id, synonym = symbol),
  by = c("taxon_id", "synonym")
)
nrow(synonyms_df)

# Remove synonyms that match the symbol of another gene in the same species,
# ignoring case
# For example, mouse FAS is a synonym of Fasn, but users mean Fas
# Keep case variants of the gene's own symbol, such as ABCD for Abcd
symbol_keys <- transmute(
  genes_df,
  taxon_id,
  key = str_to_lower(symbol),
  symbol_agr_id = agr_id
)
symbol_clashes <- synonyms_df |>
  mutate(key = str_to_lower(synonym)) |>
  inner_join(
    symbol_keys,
    by = c("taxon_id", "key"),
    relationship = "many-to-many"
  ) |>
  filter(symbol_agr_id != agr_id)
synonyms_df <- anti_join(
  synonyms_df,
  symbol_clashes,
  by = c("agr_id", "synonym")
)
nrow(synonyms_df)

# Keep synonyms that identify only one gene in each species
synonyms_df <- filter(synonyms_df, n() == 1, .by = c(taxon_id, synonym))
nrow(synonyms_df)

# Keep the relevant columns
synonyms_df <- select(synonyms_df, agr_id, synonym)
nrow(synonyms_df)

# Prepare the final package tables ----

# Keep only the columns used by the package
genes_df <- select(
  genes_df,
  taxon_id,
  symbol,
  ncbi_id,
  ensembl_id,
  agr_id
)
nrow(genes_df)

# Check that pair symbols match the gene table and remove the duplicate columns
for (side in c("g1", "g2")) {
  index <- match(orthologs_df[[paste0(side, "_agr_id")]], genes_df$agr_id)
  stopifnot(!anyNA(index))
  stopifnot(identical(
    orthologs_df[[paste0(side, "_symbol")]],
    genes_df$symbol[index]
  ))
}
orthologs_df[c("g1_symbol", "g2_symbol")] <- NULL
nrow(orthologs_df)

# Genes without a NCBI ID
genes_missing_ncbi <- filter(genes_df, is.na(ncbi_id))
genes_missing_ncbi
stopifnot(nrow(genes_missing_ncbi) < 10)
stopifnot(n_distinct(genes_missing_ncbi$taxon_id) <= 2)

# Genes with duplicated NCBI ID
genes_duplicated_ncbi <- genes_df |>
  filter(!is.na(ncbi_id)) |>
  filter(n() > 1, .by = ncbi_id) |>
  arrange(ncbi_id, symbol)
genes_duplicated_ncbi
stopifnot(nrow(genes_duplicated_ncbi) < 10)
stopifnot(n_distinct(genes_duplicated_ncbi$taxon_id) <= 1)

# Exclude problematic genes from the final tables
excluded <- c(genes_missing_ncbi$agr_id, genes_duplicated_ncbi$agr_id)
orthologs_df <- orthologs_df[
  !orthologs_df$g1_agr_id %in% excluded &
    !orthologs_df$g2_agr_id %in% excluded,
]
synonyms_df <- synonyms_df[!synonyms_df$agr_id %in% excluded, ]
genes_df <- genes_df[!genes_df$agr_id %in% excluded, ]
stopifnot(!anyNA(genes_df$ncbi_id))
stopifnot(!anyDuplicated(genes_df$ncbi_id))
stopifnot(is.integer(genes_df$ncbi_id))
stopifnot(all(genes_df$ncbi_id > 0))

# NCBI IDs are the primary gene IDs in the saved tables
# Replace AGR gene IDs with NCBI gene IDs in the orthologs table
for (side in c("g1", "g2")) {
  index <- match(orthologs_df[[paste0(side, "_agr_id")]], genes_df$agr_id)
  stopifnot(!anyNA(index))
  orthologs_df[[paste0(side, "_agr_id")]] <- genes_df$ncbi_id[index]
}
orthologs_df <- rename(
  orthologs_df,
  g1_ncbi_id = g1_agr_id,
  g2_ncbi_id = g2_agr_id
)
stopifnot(!anyDuplicated(orthologs_df[c("g1_ncbi_id", "g2_ncbi_id")]))

# Replace AGR gene IDs with NCBI gene IDs in the synonyms table
synonyms_df <- left_join(
  synonyms_df,
  select(genes_df, agr_id, ncbi_id),
  by = "agr_id"
)
synonyms_df <- select(synonyms_df, ncbi_id, synonym)
stopifnot(!anyNA(synonyms_df$ncbi_id))

genes_df <- arrange(genes_df, taxon_id, ncbi_id)
genes_df <- as.data.frame(genes_df)
rownames(genes_df) <- NULL

# A symbol must not be the NCBI ID of a different gene in the same species
id_keys <- paste(genes_df$taxon_id, genes_df$ncbi_id)
symbol_keys <- paste(genes_df$taxon_id, genes_df$symbol)
stopifnot(!any(symbol_keys %in% id_keys & symbol_keys != id_keys))

synonyms_df <- arrange(synonyms_df, ncbi_id, synonym)
synonyms_df <- as.data.frame(synonyms_df)
rownames(synonyms_df) <- NULL

orthologs_df <- select(
  orthologs_df,
  g1_taxon_id,
  g2_taxon_id,
  g1_ncbi_id,
  g2_ncbi_id,
  g1_best_match,
  g2_best_match,
  algorithms,
  algorithms_n
)
# Few algorithm combinations repeat often, so a factor compresses better
# The package functions convert it back to character
# Radix order gives the same levels in all locales
orthologs_df$algorithms <- factor(
  orthologs_df$algorithms,
  levels = sort(unique(orthologs_df$algorithms), method = "radix")
)
# Keep the pairs of each gene together because this compresses better
orthologs_df <- arrange(
  orthologs_df,
  g1_taxon_id,
  g1_ncbi_id,
  g2_taxon_id,
  g2_ncbi_id
)
orthologs_df <- as.data.frame(orthologs_df)
rownames(orthologs_df) <- NULL

species_df <- arrange(species_df, taxon_id, name, name_class)
species_df <- as.data.frame(species_df)
rownames(species_df) <- NULL

# Check that the final tables use the same genes and species
ortholog_ids <- c(orthologs_df$g1_ncbi_id, orthologs_df$g2_ncbi_id)
stopifnot(setequal(genes_df$ncbi_id, ortholog_ids))
stopifnot(all(synonyms_df$ncbi_id %in% genes_df$ncbi_id))
for (side in c("g1", "g2")) {
  index <- match(orthologs_df[[paste0(side, "_ncbi_id")]], genes_df$ncbi_id)
  stopifnot(identical(
    genes_df$taxon_id[index],
    orthologs_df[[paste0(side, "_taxon_id")]]
  ))
}
stopifnot(setequal(species_df$taxon_id, genes_df$taxon_id))

# Save the package tables ----

usethis::use_data(
  genes_df,
  synonyms_df,
  orthologs_df,
  species_df,
  internal = TRUE,
  overwrite = TRUE,
  compress = "xz"
)

# Delete downloaded source files after saving the package tables
unlink(source_dir, recursive = TRUE)
