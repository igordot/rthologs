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
