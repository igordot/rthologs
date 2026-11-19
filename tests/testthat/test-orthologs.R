test_that("symbols, NCBI IDs, and Ensembl IDs give the same pairs", {
  out <- orthologs("TP53", source = "human", target = "mouse")
  expect_identical(out$input, "TP53")
  expect_identical(out$target_symbol, "Trp53")
  expect_identical(out$source_ensembl_id, "ENSG00000141510")
  expect_identical(out$target_ensembl_id, "ENSMUSG00000059552")
  for (input in c("7157", "ENSG00000141510")) {
    by_id <- orthologs(input, "human", "mouse")
    expect_identical(by_id$input, input)
    expect_identical(by_id[-1], out[-1])
  }
  # Each different input gets one row
  all <- orthologs(
    c("TP53", "7157", "ENSG00000141510", "TP53"),
    "human",
    "mouse"
  )
  expect_identical(all$input, c("TP53", "7157", "ENSG00000141510"))
  expect_identical(unique(all[-1]), out[-1])
  reverse <- orthologs("Trp53", "mouse", "human")
  expect_identical(reverse$target_ncbi_id, out$source_ncbi_id)
  expect_identical(reverse$source_ncbi_id, out$target_ncbi_id)
  expect_identical(orthologs("Trp53", "mouse", "rat")$target_symbol, "Tp53")
})

test_that("taxon IDs give the same pairs as species names", {
  expect_identical(
    orthologs("TP53", 9606, "10090"),
    orthologs("TP53", "human", "mouse")
  )
})

test_that("synonyms give the same pairs as the current symbol", {
  out <- orthologs("TP53", "human", "mouse")
  by_synonym <- orthologs("LFS1", "human", "mouse")
  expect_identical(by_synonym$input, "LFS1")
  expect_identical(by_synonym$source_symbol, "TP53")
  expect_identical(by_synonym[-1], out[-1])
})

test_that("synonyms = FALSE keeps only direct matches", {
  expect_no_warning(
    out <- orthologs(c("PTPRC", "LFS1"), "human", "mouse", synonyms = FALSE)
  )
  expect_identical(out$source_symbol, "PTPRC")
})

test_that("synonyms of a different species do not match", {
  expect_identical(orthologs("Ly-5", "mouse", "human")$target_symbol, "PTPRC")
  expect_warning(out <- orthologs("Ly-5", "human", "mouse"))
  expect_identical(nrow(out), 0L)
})
