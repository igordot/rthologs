test_that("species names resolve to the same taxon", {
  expect_identical(species("mouse"), species(" Mus musculus "))
  expect_identical(unique(species("MOUSE")$taxon_id), 10090L)
})

test_that("taxon IDs resolve to the same taxon as names", {
  expect_identical(species(10090), species("mouse"))
  expect_identical(species(10090L), species("mouse"))
  expect_identical(species(" 10090 "), species("mouse"))
})

test_that("invalid species names cause errors", {
  expect_error(species("not a species"))
  expect_error(species(c("mouse", "rat")))
  expect_error(species(NA_character_))
  expect_error(species(99999))
  expect_error(species(10090.5))
  expect_error(species(TRUE))
})

test_that("species() returns the complete catalog", {
  catalog <- species()
  expect_named(catalog, c("taxon_id", "name", "name_class"))
  taxa <- c(6239L, 7227L, 7955L, 8364L, 9606L, 10090L, 10116L, 559292L)
  expect_setequal(unique(catalog$taxon_id), taxa)
  for (id in taxa) {
    listed <- catalog[catalog$taxon_id == id, ]
    expect_gte(nrow(listed), 2L)
    expect_true(dplyr::setequal(listed, species(id)), info = as.character(id))
  }
})

test_that("all name records resolve to one species", {
  records <- unique(species()[c("taxon_id", "name")])
  for (i in seq_len(nrow(records))) {
    expect_identical(
      unique(species(records$name[i])$taxon_id),
      records$taxon_id[i]
    )
  }
})
