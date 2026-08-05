# Tests for CreateATACObjects()/CreateATACObjectsFilter()'s argument
# validation only -- `genome`, `object_names`, and (pre-existing)
# add_treatment/treatment consistency, and (for the Filter variant)
# filter/interactive consistency. All of these are checked before any file
# is read or any package requireNamespace()'d, so they run with bogus
# `data_dirs` and no ATAC-specific packages installed. The read+build
# pipeline itself is not covered here -- see the note in
# test-object-construction.R for why.

test_that("CreateATACObjects errors on an invalid `genome`", {
  expect_error(
    CreateATACObjects(data_dirs = "nonexistent_dir", genome = "hg19"),
    "arg"
  )
})

test_that("CreateATACObjects errors when treatment is set but add_treatment is FALSE", {
  expect_error(
    CreateATACObjects(data_dirs = "nonexistent_dir", treatment = "A"),
    "add_treatment"
  )
})

test_that("CreateATACObjects errors when object_names doesn't match length(data_dirs)", {
  expect_error(
    CreateATACObjects(data_dirs = c("dir1", "dir2"), object_names = "only_one"),
    "object_names"
  )
})

test_that("CreateATACObjectsFilter errors on an invalid `genome`", {
  expect_error(
    CreateATACObjectsFilter(data_dirs = "nonexistent_dir", genome = "hg19"),
    "arg"
  )
})

test_that("CreateATACObjectsFilter errors when interactive = TRUE but filter = FALSE", {
  expect_error(
    CreateATACObjectsFilter(data_dirs = "nonexistent_dir", filter = FALSE, interactive = TRUE)
  )
})

test_that("CreateATACObjectsFilter errors when object_names doesn't match length(data_dirs)", {
  expect_error(
    CreateATACObjectsFilter(data_dirs = c("dir1", "dir2"), object_names = "only_one"),
    "object_names"
  )
})
