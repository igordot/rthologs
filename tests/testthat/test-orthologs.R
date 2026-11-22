test_that("symbols, NCBI IDs, and Ensembl IDs give the same pairs", {
  out <- orthologs("TP53", source = "human", target = "mouse", details = TRUE)
  expect_identical(out$input, "TP53")
  expect_identical(out$target_symbol, "Trp53")
  expect_identical(out$source_ensembl_id, "ENSG00000141510")
  expect_identical(out$target_ensembl_id, "ENSMUSG00000059552")
  for (input in c("7157", "ENSG00000141510")) {
    by_id <- orthologs(input, "human", "mouse", details = TRUE)
    expect_identical(by_id$input, input)
    expect_identical(by_id[-1], out[-1])
  }
  # Numbers are searched as text
  by_id <- orthologs("7157", "human", "mouse")
  expect_identical(orthologs(7157L, "human", "mouse"), by_id)
  expect_identical(orthologs(7157, "human", "mouse"), by_id)
  # Each different input gets one row
  all <- orthologs(
    c("TP53", "7157", "ENSG00000141510", "TP53"),
    "human",
    "mouse",
    details = TRUE
  )
  expect_identical(all$input, c("TP53", "7157", "ENSG00000141510"))
  expect_identical(unique(all[-1]), out[-1])
  reverse <- orthologs("Trp53", "mouse", "human", details = TRUE)
  expect_identical(reverse$target_ncbi_id, out$source_ncbi_id)
  expect_identical(reverse$source_ncbi_id, out$target_ncbi_id)
  expect_identical(orthologs("Trp53", "mouse", "rat")$target_symbol, "Tp53")
})

test_that("numeric NCBI IDs match when R prints them in scientific notation", {
  expected <- orthologs("500000", "rat", "human")
  expect_identical(expected$source_symbol, "Cldn12")
  expect_identical(expected$target_symbol, "CLDN12")
  expect_identical(orthologs(500000, "rat", "human"), expected)
  expect_identical(orthologs(500000L, "rat", "human"), expected)
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
    out <- orthologs(
      c("PTPRC", "LFS1"),
      "human",
      "mouse",
      synonyms = FALSE
    )
  )
  expect_identical(out$source_symbol, "PTPRC")
})

test_that("synonyms of a different species do not match", {
  expect_identical(orthologs("Ly-5", "mouse", "human")$target_symbol, "PTPRC")
  expect_warning(out <- orthologs("Ly-5", "human", "mouse"))
  expect_identical(nrow(out), 0L)
})

test_that("synonyms do not match other symbols in another case", {
  # Each input is a synonym of a different gene than the one with that symbol
  mm <- c(
    "FAS",
    "EPO",
    "SLC7A1",
    "BCR",
    "TEP1",
    "TPO",
    "PGF",
    "HP",
    "GRIN1",
    "DCC"
  )
  expect_warning(
    out <- orthologs(mm, "mouse", "human")
  )
  expect_identical(nrow(out), 0L)
  expect_warning(
    out <- orthologs("Fli1", "human", "mouse")
  )
  expect_identical(nrow(out), 0L)
})

# Expected gene pairs
test_that("the five marker genes return their mouse orthologs", {
  out <- orthologs(
    c("PTPRC", "EPCAM", "CD34", "CD3D", "CD19"),
    "human",
    "mouse"
  )
  expect_identical(
    out[c("source_symbol", "target_ncbi_id", "target_symbol")],
    data.frame(
      source_symbol = c("PTPRC", "EPCAM", "CD34", "CD3D", "CD19"),
      target_ncbi_id = c("19264", "17075", "12490", "12500", "12478"),
      target_symbol = c("Ptprc", "Epcam", "Cd34", "Cd3d", "Cd19")
    )
  )
})

test_that("PTPRC has the expected rat and zebrafish genes", {
  rn <- orthologs("PTPRC", "human", "rat")
  expect_identical(rn$target_ncbi_id, "24699")
  expect_identical(rn$target_symbol, "Ptprc")
  dr <- orthologs("PTPRC", "human", "zebrafish")
  expect_identical(dr$target_ncbi_id, "559154")
  expect_identical(dr$target_symbol, "ptprc")
  # A pair without a human gene
  expect_identical(orthologs("Ptprc", "mouse", "rat")$target_ncbi_id, "24699")
  expect_identical(orthologs("Ptprc", "rat", "mouse")$target_ncbi_id, "19264")
})

test_that("known difficult genes give the expected pairs", {
  expect_no_warning(
    ensa <- orthologs(
      c("ENSG00000083622", "ENSG00000145063", "ENSA"),
      "human",
      "mouse"
    )
  )
  expect_identical(ensa$input, "ENSA")
  expect_identical(ensa$target_symbol, "Ensa")
  ndst2 <- orthologs("NDST2", "human", "mouse")
  expect_identical(ndst2$target_symbol, "Ndst2")
  expect_identical(
    ndst2[-1],
    orthologs("ENSG00000166507", "human", "mouse")[-1]
  )
  # LILRB3 has no mouse pair
  expect_warning(
    lilrb3 <- orthologs("LILRB3", "human", "mouse")
  )
  expect_identical(nrow(lilrb3), 0L)
  lilrb3 <- orthologs("LILRB3", "human", "rat", details = TRUE)
  expect_identical(lilrb3$target_symbol, "Pirb")
  expect_identical(lilrb3$support_n, 1L)
  # Rat Pals1: the newer Ensembl files give this Ensembl ID the name MPP5
  pals1 <- orthologs("ENSRNOG00000008788", "rat", "human")
  expect_identical(pals1$source_symbol, "Pals1")
  expect_identical(pals1$target_symbol, "PALS1")
  # Fly 14-3-3 epsilon: the symbol has the Greek letter, written as ε
  # Some sources write 14-3-3epsilon, which is a synonym
  dm <- "Drosophila melanogaster"
  greek <- orthologs("14-3-3ε", dm, "human", details = TRUE)
  expect_identical(greek$source_ncbi_id, "42186")
  expect_identical(greek$target_symbol, "YWHAE")
  epsilon <- orthologs("14-3-3epsilon", dm, "human", details = TRUE)
  expect_identical(epsilon[-1], greek[-1])
  # TAS2R43 and TAS2R45 have the same Ensembl ID
  # AGR gives it for TAS2R45 and HGNC gives it for TAS2R43
  shared <- orthologs("ENSG00000255374", "human", "mouse")
  expect_identical(shared$source_symbol, c("TAS2R43", "TAS2R45"))
  expect_identical(shared$target_symbol, rep("Tas2r120", 2))
})

test_that("BGLAP keeps three equal best mouse matches", {
  out <- orthologs("BGLAP", "human", "mouse", details = TRUE)
  expect_identical(out$target_ncbi_id, c("12095", "12096", "12097"))
  expect_identical(out$target_symbol, c("Bglap3", "Bglap", "Bglap2"))
  expect_identical(out$support_n, rep(10L, 3))
  expect_identical(
    out,
    orthologs("BGLAP", "human", "mouse", best = FALSE, details = TRUE)
  )
})

test_that("SYCP3 has one mouse match by default and 76 without best", {
  expect_identical(orthologs("SYCP3", "human", "mouse")$target_symbol, "Sycp3")
  expect_identical(
    nrow(orthologs("SYCP3", "human", "mouse", best = FALSE)),
    76L
  )
  expect_identical(orthologs("Sycp3", "mouse", "human")$target_symbol, "SYCP3")
})

test_that("GAPDH filters remove mouse matches with lower support", {
  out <- orthologs("GAPDH", "human", "mouse", best = FALSE, details = TRUE)
  expect_identical(out$target_ncbi_id, c("100039229", "100042427", "14433"))
  expect_identical(out$target_symbol, c("Gapdhrt2", "Gapdhrt", "Gapdh"))
  expect_identical(out$support_n, c(9L, 9L, 10L))
  best <- orthologs("GAPDH", "human", "mouse")
  expect_identical(best$target_ncbi_id, "14433")
  expect_identical(
    best,
    orthologs("GAPDH", "human", "mouse", min_support = 10, best = FALSE)
  )
  expect_identical(
    out,
    orthologs(
      "GAPDH",
      "human",
      "mouse",
      min_support = 9,
      best = FALSE,
      details = TRUE
    )
  )
})

test_that("`best` uses the stored best-match flags, not support", {
  pairs <- data.frame(
    g1_taxon_id = 1L,
    g1_ncbi_id = c("A", "A", "A", "B", "C"),
    g1_symbol = c("same", "same", "same", "same", "123"),
    g2_taxon_id = 2L,
    g2_ncbi_id = c("X", "Y", "Z", "Z", "W"),
    g2_symbol = c("x", "y", "z", "z", "w"),
    algorithms = c("a|b|c", "a|b|c", "a|b", "a", "a"),
    algorithms_n = c(3L, 3L, 2L, 1L, 1L),
    g1_best_match = c(TRUE, FALSE, TRUE, TRUE, FALSE),
    g2_best_match = c(TRUE, TRUE, TRUE, FALSE, TRUE)
  )
  gene_table <- unique(rbind(
    setNames(pairs[c("g1_ncbi_id", "g1_symbol")], c("ncbi_id", "symbol")),
    setNames(pairs[c("g2_ncbi_id", "g2_symbol")], c("ncbi_id", "symbol"))
  ))
  gene_table$ensembl_id <- NA_character_
  pairs[c("g1_symbol", "g2_symbol")] <- NULL
  forward <- lookup_orthologs("same", 1L, 2L, 1, TRUE, pairs, gene_table)
  expect_identical(forward$target_ncbi_id, c("X", "Z", "Z"))
  expect_identical(forward$source_ncbi_id, c("A", "A", "B"))
  reverse <- lookup_orthologs("Z", 2L, 1L, 1, TRUE, pairs, gene_table)
  expect_identical(reverse$target_ncbi_id, "A")
  expect_identical(
    nrow(lookup_orthologs("same", 1L, 2L, 1, FALSE, pairs, gene_table)),
    4L
  )
  expect_identical(
    lookup_orthologs("123", 1L, 2L, 1, FALSE, pairs, gene_table)$target_ncbi_id,
    "W"
  )
  expect_identical(
    lookup_orthologs("A", 1L, 2L, 3, FALSE, pairs, gene_table)$target_ncbi_id,
    c("X", "Y")
  )
})

# Build the expected orthologs() result from the saved pairs of one pair
expected_pairs <- function(pairs, s, t) {
  gene_columns <- function(side, prefix) {
    ids <- pairs[[paste0(side, "_ncbi_id")]]
    rows <- match(ids, genes_df$ncbi_id)
    out <- data.frame(
      taxon_id = pairs[[paste0(side, "_taxon_id")]],
      symbol = genes_df$symbol[rows],
      ncbi_id = as.character(ids),
      ensembl_id = genes_df$ensembl_id[rows]
    )
    names(out) <- paste0(prefix, "_", names(out))
    out
  }
  out <- cbind(
    gene_columns(s, "source"),
    gene_columns(t, "target"),
    support = as.character(pairs$algorithms),
    support_n = pairs$algorithms_n,
    best = pairs[[paste0(s, "_best_match")]]
  )
  out <- out[order(out$source_ncbi_id, out$target_ncbi_id), ]
  rownames(out) <- NULL
  out
}

test_that("`best` gives the stored pairs for every species pair", {
  combinations <- unique(orthologs_df[c("g1_taxon_id", "g2_taxon_id")])
  for (i in seq_len(nrow(combinations))) {
    pairs <- orthologs_df[
      orthologs_df$g1_taxon_id == combinations$g1_taxon_id[i] &
        orthologs_df$g2_taxon_id == combinations$g2_taxon_id[i],
    ]
    for (side in list(c("g1", "g2"), c("g2", "g1"))) {
      expected <- expected_pairs(pairs, side[1], side[2])
      best <- expected$best
      expected$best <- NULL
      inputs <- unique(expected$source_ncbi_id)
      source <- expected$source_taxon_id[1]
      target <- expected$target_taxon_id[1]
      out <- orthologs(inputs, source, target, best = FALSE, details = TRUE)
      # The inputs are NCBI IDs, so `input` repeats the source ID
      expect_identical(out$input, out$source_ncbi_id)
      expect_identical(out[-1], expected)
      expected_best <- expected[best, ]
      rownames(expected_best) <- NULL
      out <- orthologs(inputs, source, target, details = TRUE)
      expect_identical(out[-1], expected_best)
    }
  }
})

test_that("`orthologs_only = FALSE` keeps inputs without orthologs", {
  out <- orthologs(c("NOTAGENE", "TP53"), "human", "mouse")
  kept <- orthologs(
    c("NOTAGENE", "TP53"),
    "human",
    "mouse",
    orthologs_only = FALSE
  )
  expect_identical(kept$input, c("NOTAGENE", "TP53"))
  expect_true(all(is.na(kept[1, -1])))
  expect_identical(kept[2, -1], out[1, -1], ignore_attr = TRUE)
  expect_identical(
    vapply(kept, typeof, character(1)),
    vapply(out, typeof, character(1))
  )
  # Inputs that the filters removed also get a row
  expect_warning(
    kept <- orthologs(
      "TP53",
      "human",
      "mouse",
      min_support = 100,
      orthologs_only = FALSE
    )
  )
  expect_identical(kept$input, "TP53")
  expect_true(is.na(kept$target_symbol))
})

test_that("only missing values and empty results cause warnings", {
  expect_warning(out <- orthologs(c("TP53", NA_character_), "human", "mouse"))
  expect_identical(out$target_symbol, "Trp53")
  # Filters that remove all pairs of a synonym also cause a warning
  expect_warning(out <- orthologs("LFS1", "human", "mouse", min_support = 100))
  expect_identical(nrow(out), 0L)
  expect_warning(
    out <- orthologs(
      c("TP53", "not_a_gene"),
      "human",
      "mouse",
      min_support = 100
    )
  )
  expect_identical(nrow(out), 0L)
  expect_no_warning(empty <- orthologs(character(), "human", "mouse"))
  expect_identical(out, empty)
  expect_identical(
    names(empty),
    c(
      "input",
      "source_symbol",
      "target_symbol",
      "target_ncbi_id",
      "target_ensembl_id"
    )
  )
  empty <- orthologs(character(), "human", "mouse", details = TRUE)
  expect_identical(
    vapply(empty, typeof, character(1)),
    c(
      input = "character",
      source_taxon_id = "integer",
      source_symbol = "character",
      source_ncbi_id = "character",
      source_ensembl_id = "character",
      target_taxon_id = "integer",
      target_symbol = "character",
      target_ncbi_id = "character",
      target_ensembl_id = "character",
      support = "character",
      support_n = "integer"
    )
  )
})

test_that("invalid input and unknown species cause errors", {
  expect_error(orthologs(list("TP53"), "human", "mouse"))
  expect_error(orthologs(matrix("TP53"), "human", "mouse"))
  expect_error(orthologs("TP53", "human", "mouse", best = NA))
  expect_error(orthologs("TP53", "human", "mouse", best = c(TRUE, FALSE)))
  expect_error(orthologs("TP53", "human", "mouse", min_support = 1.5))
  expect_error(orthologs("TP53", "human", "mouse", min_support = 0))
  expect_error(orthologs("TP53", "human", "mouse", synonyms = NA))
  expect_error(orthologs("TP53", "human", "mouse", orthologs_only = NA))
  expect_error(orthologs("TP53", "human", "mouse", details = NA))
  expect_error(orthologs("TP53", NULL, "mouse"))
  expect_error(orthologs("TP53", "human", "human"))
  expect_error(orthologs("Trp53", "mouse", "Xenopus laevis"))
})

# Get each pair in both search directions with the species names
read_searches <- function() {
  cases <- read.delim(test_path("fixtures", "difficult-pairs.tsv"))
  reverse <- cases[c(6:10, 1:5)]
  names(reverse) <- names(cases)
  cases <- rbind(cases, reverse)
  species_names <- species_df[species_df$name_class == "species name", ]
  cases$source <- species_names$name[match(
    cases$taxon_id,
    species_names$taxon_id
  )]
  cases$target <- species_names$name[match(
    cases$partner_taxon_id,
    species_names$taxon_id
  )]
  cases
}

# Search each gene by one input type and return the genes with no results
search_failures <- function(column) {
  cases <- read_searches()
  cases <- cases[!is.na(cases[[column]]), ]
  found <- vapply(
    seq_len(nrow(cases)),
    function(i) {
      out <- orthologs(
        cases[[column]][i],
        cases$source[i],
        cases$target[i],
        details = TRUE
      )
      nrow(out) > 0L && all(out$source_ncbi_id == cases$ncbi_id[i])
    },
    logical(1)
  )
  cases$symbol[!found]
}

test_that("a search by symbol finds the expected gene", {
  expect_identical(search_failures("symbol"), character())
})

test_that("a search by NCBI ID finds the expected gene", {
  expect_identical(search_failures("ncbi_id"), character())
})

test_that("a search by Ensembl ID finds the expected gene", {
  expect_identical(search_failures("ensembl_id"), character())
})

test_that("each difficult pair is present in both search directions", {
  cases <- read_searches()
  for (i in seq_len(nrow(cases))) {
    out <- orthologs(
      cases$ncbi_id[i],
      cases$source[i],
      cases$target[i],
      best = FALSE,
      details = TRUE
    )
    targets <- out$target_ncbi_id[out$source_ncbi_id == cases$ncbi_id[i]]
    expect_in(as.character(cases$partner_ncbi_id[i]), targets)
  }
})

test_that("MSigDB M45000 human genes match MM1334", {
  # `source_id` is a human Ensembl ID
  # MM1334 is the mouse version of M45000
  hs_set <- read.delim(test_path("fixtures", "genes-human-human-M45000.tsv"))
  mm_set <- read.delim(test_path("fixtures", "genes-mouse-human-MM1334.tsv"))
  in_mm <- hs_set$source_id %in% mm_set$source_id
  out <- orthologs(hs_set$symbol[in_mm], "human", "mouse", details = TRUE)
  expect_lte(sum(!hs_set$ncbi_id[in_mm] %in% out$source_ncbi_id), 1)
  expect_lte(sum(!hs_set$source_id[in_mm] %in% out$source_ensembl_id), 1)
  # The genes that MSigDB does not map to a mouse gene have no orthologs
  out <- orthologs(hs_set$symbol[!in_mm], "human", "mouse")
  expect_lt(length(unique(out$input)) / sum(!in_mm), 0.05)
})

test_that("MSigDB MM17506 finds the same genes for human and rat", {
  # `source_id` is a mouse Ensembl ID
  gene_set <- read.delim(test_path("fixtures", "genes-mouse-mouse-MM17506.tsv"))
  hs <- orthologs(gene_set$symbol, "mouse", "human", details = TRUE)
  rn <- orthologs(gene_set$symbol, "mouse", "rat", details = TRUE)
  expect_lte(sum(!gene_set$source_id %in% hs$source_ensembl_id), 1)
  hs <- intersect(gene_set$ncbi_id, hs$source_ncbi_id)
  rn <- intersect(gene_set$ncbi_id, rn$source_ncbi_id)
  expect_gt(length(hs) / length(unique(gene_set$ncbi_id)), 0.95)
  expect_gt(length(rn) / length(unique(gene_set$ncbi_id)), 0.95)
  expect_gt(length(intersect(hs, rn)) / length(union(hs, rn)), 0.95)
})

test_that("MSigDB MM3637 finds the same genes for three species", {
  gene_set <- read.delim(test_path("fixtures", "genes-mouse-mouse-MM3637.tsv"))
  hs <- orthologs(gene_set$symbol, "mouse", "human", details = TRUE)
  rn <- orthologs(gene_set$symbol, "mouse", "rat", details = TRUE)
  dr <- orthologs(gene_set$symbol, "mouse", "zebrafish", details = TRUE)
  hs <- intersect(gene_set$ncbi_id, hs$source_ncbi_id)
  rn <- intersect(gene_set$ncbi_id, rn$source_ncbi_id)
  dr <- intersect(gene_set$ncbi_id, dr$source_ncbi_id)
  expect_gt(length(hs) / length(unique(gene_set$ncbi_id)), 0.8)
  expect_gt(length(rn) / length(unique(gene_set$ncbi_id)), 0.8)
  expect_gt(length(dr) / length(unique(gene_set$ncbi_id)), 0.7)
  expect_gt(length(intersect(hs, rn)) / length(union(hs, rn)), 0.95)
  all <- intersect(intersect(hs, rn), dr)
  any <- union(union(hs, rn), dr)
  expect_gt(length(all) / length(any), 0.7)
})

test_that("MSigDB MM1334 human Ensembl IDs give the mouse NCBI IDs", {
  # `source_id` is a human Ensembl ID
  gene_set <- read.delim(test_path("fixtures", "genes-mouse-human-MM1334.tsv"))
  out <- orthologs(gene_set$source_id, "human", "mouse")
  same <- paste(gene_set$source_id, gene_set$ncbi_id) %in%
    paste(out$input, out$target_ncbi_id)
  expect_lte(sum(!same), 1)
  # The reverse conversion gives the MSigDB human Ensembl IDs for most genes
  reverse <- orthologs(gene_set$symbol, "mouse", "human")
  same <- paste(gene_set$symbol, gene_set$source_id) %in%
    paste(reverse$input, reverse$target_ensembl_id)
  expect_gt(mean(same), 0.8)
})

test_that("MSigDB M46418 rat symbols give the human NCBI IDs", {
  # `source_id` is a rat symbol
  gene_set <- read.delim(test_path("fixtures", "genes-human-rat-M46418.tsv"))
  out <- orthologs(gene_set$source_id, "rat", "human")
  same <- paste(gene_set$source_id, gene_set$ncbi_id) %in%
    paste(out$input, out$target_ncbi_id)
  expect_gt(mean(same), 0.95)
})

test_that("MSigDB MM16096 rat symbols give the mouse NCBI IDs", {
  # `source_id` is a rat symbol
  gene_set <- read.delim(test_path("fixtures", "genes-mouse-rat-MM16096.tsv"))
  out <- orthologs(gene_set$source_id, "rat", "mouse")
  same <- paste(gene_set$source_id, gene_set$ncbi_id) %in%
    paste(out$input, out$target_ncbi_id)
  expect_gt(mean(same), 0.95)
})
