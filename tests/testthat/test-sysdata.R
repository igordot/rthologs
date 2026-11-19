test_that("supported species exclude X. laevis and have every species pair", {
  taxa <- sort(unique(species_df$taxon_id))
  expect_identical(
    taxa,
    c(6239L, 7227L, 7955L, 8364L, 9606L, 10090L, 10116L, 559292L)
  )
  combinations <- unique(orthologs_df[c("g1_taxon_id", "g2_taxon_id")])
  expect_identical(nrow(combinations), as.integer(choose(length(taxa), 2)))
  partners <- table(c(combinations$g1_taxon_id, combinations$g2_taxon_id))
  expect_true(all(partners == length(taxa) - 1L))
  expect_false(8355L %in% genes_df$taxon_id)
})

test_that("internal best-score flags are logical without missing values", {
  expect_type(orthologs_df$g1_best_match, "logical")
  expect_type(orthologs_df$g2_best_match, "logical")
  expect_identical(
    anyNA(orthologs_df[c("g1_best_match", "g2_best_match")]),
    FALSE
  )
})

test_that("the gene table covers each ortholog gene once", {
  ids <- unique(c(orthologs_df$g1_ncbi_id, orthologs_df$g2_ncbi_id))
  expect_named(
    genes_df,
    c("taxon_id", "symbol", "ncbi_id", "ensembl_id", "agr_id")
  )
  expect_type(ids, "integer")
  expect_type(genes_df$ncbi_id, "integer")
  expect_setequal(genes_df$ncbi_id, ids)
  expect_identical(anyDuplicated(genes_df$agr_id), 0L)
  expect_identical(anyDuplicated(genes_df$ncbi_id), 0L)
  expect_false(anyNA(genes_df[c("taxon_id", "agr_id", "ncbi_id", "symbol")]))
  expect_true(all(genes_df$ncbi_id > 0L))
  expect_length(
    intersect(genes_df$ncbi_id, c(553461L, 791178L, 569678L, 565755L)),
    0L
  )
})

test_that("Ensembl IDs are single values or missing", {
  id <- genes_df$ensembl_id
  expect_type(id, "character")
  expect_false(any(grepl("|", id, fixed = TRUE)))
  expect_true(all(nzchar(id[!is.na(id)])))
  expect_lt(mean(is.na(id)), 0.04)
})

test_that("NCBI and Ensembl fill missing identifiers", {
  expect_true(all(c(553997L, 100150068L, 101884914L) %in% genes_df$ncbi_id))
  gene <- genes_df[genes_df$ncbi_id == 101884914L, ]
  expect_identical(gene$ensembl_id, "ENSDARG00000097875")
  gene <- genes_df[genes_df$ncbi_id == 100150068L, ]
  expect_identical(gene$ensembl_id, "ENSDARG00180002076")
  gene <- genes_df[genes_df$ncbi_id == 853033L, ]
  expect_identical(gene$ensembl_id, "YGR132C")
})

test_that("genes with multiple AGR NCBI IDs get the ID with the matching symbol", {
  expected <- c(
    "SGD:S000005659" = 854301L,
    "SGD:S000002793" = 851993L,
    "MGI:2444426" = 210027L,
    "SGD:S000007274" = 854584L
  )
  rows <- match(names(expected), genes_df$agr_id)
  expect_identical(genes_df$ncbi_id[rows], unname(expected))
})

test_that("HGNC fills human genes without a single AGR Ensembl ID", {
  # AAAS, AADACL2, ABCC1, CMAHP, AQP7P3
  ids <- c("8086", "344752", "4363", "8418", "441432")
  ens <- genes_df[match(ids, genes_df$ncbi_id), ]
  expect_identical(
    ens$ensembl_id,
    c(
      "ENSG00000094914",
      "ENSG00000197953",
      "ENSG00000103222",
      "ENSG00000168405",
      "ENSG00000279688"
    )
  )
})

# Bioconductor OrgDb gives another Ensembl ID for these human genes
test_that("human genes with unconfirmed Ensembl IDs keep the current IDs", {
  # FOLH1B, GABARAPL3, KLRC4-KLRK1, OR4Q2
  ids <- c(219595L, 23766L, 100528032L, 390432L)
  genes <- genes_df[match(ids, genes_df$ncbi_id), ]
  expect_identical(
    genes$ensembl_id,
    c(
      "ENSG00000134612",
      "ENSG00000238244",
      "ENSG00000255819",
      "ENSG00000196383"
    )
  )
})

test_that("Ensembl and NCBI select among multiple AGR Ensembl IDs", {
  # Rat Rps10-ps22: AGR lists two IDs, and only one is in the Ensembl GTF file
  ens <- genes_df[genes_df$ncbi_id == "299823", ]
  expect_identical(ens$ensembl_id, "ENSRNOG00000049099")
  # Mouse Slc60a2-ps2: both AGR IDs are in the GTF file, and NCBI lists one
  ens <- genes_df[genes_df$ncbi_id == "100041085", ]
  expect_identical(ens$ensembl_id, "ENSMUSG00000121339")
})

test_that("synonyms do not duplicate symbols or Ensembl IDs in a species", {
  taxa <- genes_df$taxon_id[match(
    synonyms_df$ncbi_id,
    genes_df$ncbi_id
  )]
  synonym_keys <- paste(taxa, synonyms_df$synonym, sep = ":")
  symbol_keys <- paste(genes_df$taxon_id, genes_df$symbol, sep = ":")
  ensembl <- genes_df[!is.na(genes_df$ensembl_id), ]
  ensembl_keys <- paste(ensembl$taxon_id, ensembl$ensembl_id, sep = ":")
  expect_length(intersect(synonym_keys, symbol_keys), 0L)
  expect_length(intersect(synonym_keys, ensembl_keys), 0L)
})

test_that("synonyms identify at most one gene per species", {
  taxa <- genes_df$taxon_id[match(
    synonyms_df$ncbi_id,
    genes_df$ncbi_id
  )]
  synonym_keys <- paste(taxa, synonyms_df$synonym, sep = ":")
  expect_identical(anyDuplicated(synonym_keys), 0L)
})
