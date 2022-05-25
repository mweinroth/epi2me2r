library(epi2me2r)

example.amr.dir <- system.file(package = "epi2me2r", "extdata", "example_amr_data")
example.wimp.dir <- system.file(package = "epi2me2r", "extdata", "example_wimp_data")

test_that("read_in_amr_files returns a data.frame", {
  df.test <- read_in_amr_files(path.to.amr.files = example.amr.dir)
  expect_s3_class(df.test, "data.frame")
})

test_that("generate_amr_taxonomy returns a data.frame", {
  amr.count.table <- read_in_amr_files(path.to.amr.files = example.amr.dir)
  df.test <- generate_amr_taxonomy(amr.count.table)
  expect_s3_class(df.test, "data.frame")
})

test_that("amr_raw_to_phyloseq returns a phyloseq object", {
  metadata <- read.csv(system.file(package = "epi2me2r", "extdata", "example_metadata.csv"))
  df.test <- amr_raw_to_phyloseq(path.to.amr.files = example.amr.dir,
                                 metadata = metadata)
  expect_s4_class(df.test, "phyloseq")
})

test_that("amr_raw_to_metagenomeseq returns a MRexperiment object", {
  metadata <- read.csv(system.file(package = "epi2me2r", "extdata", "example_metadata.csv"))
  df.test <- amr_raw_to_metagenomeseq(path.to.amr.files = example.amr.dir,
                                      metadata = metadata)
  expect_s4_class(df.test, "MRexperiment")
})

test_that("amr_read_taxonomy returns a data.frame", {
  df.test <- amr_read_taxonomy(path.to.wimp.files = example.wimp.dir, path.to.amr.files = example.amr.dir)
  expect_s3_class(df.test, "data.frame")
})



