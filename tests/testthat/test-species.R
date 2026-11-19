test_that("species names resolve to the same taxon", {
  expect_identical(species("mouse"), species(" Mus musculus "))
  expect_identical(unique(species("MOUSE")$taxon_id), 10090L)
})

test_that("taxon IDs resolve to the same taxon as names", {
  expect_identical(species(10090), species("mouse"))
  expect_identical(species(10090L), species("mouse"))
  expect_identical(species(" 10090 "), species("mouse"))
  for (id in unique(species_df$taxon_id)) {
    expect_identical(unique(species(id)$taxon_id), id)
  }
})

test_that("invalid species names cause errors", {
  expect_error(species("not a species"))
  expect_error(species(c("mouse", "rat")))
  expect_error(species(NA_character_))
  expect_error(species(99999))
  expect_error(species(10090.5))
  expect_error(species(TRUE))
})

test_that("all name records resolve to one species", {
  expect_named(species(), c("taxon_id", "name", "name_class"))
  records <- unique(species_df[c("taxon_id", "name")])
  for (i in seq_len(nrow(records))) {
    expect_identical(
      unique(species(records$name[i])$taxon_id),
      records$taxon_id[i]
    )
  }
})
