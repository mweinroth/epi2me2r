library(epi2me2r)

example.amr.dir <- system.file(package = "epi2me2r", "extdata", "example_amr_data")
example.wimp.dir <- system.file(package = "epi2me2r", "extdata", "example_wimp_data")

test_that("read_in_wimp_files returns a data.frame", {
  df.test <- read_in_wimp_files(path.to.wimp.files = example.wimp.dir)
  expect_s3_class(df.test, "data.frame")
})

test_that("generate_wimp_taxonomy returns a data.frame", {
  wimp.count.table <- read_in_wimp_files(path.to.wimp.files = example.wimp.dir)
  df.test <- generate_wimp_taxonomy(wimp.count.table)
  expect_s3_class(df.test, "data.frame")
})

test_that("wimp_raw_to_phyloseq returns a phyloseq object", {
  metadata <- read.csv(system.file(package = "epi2me2r", "extdata", "example_metadata.csv"))
  df.test <- wimp_raw_to_phyloseq(path.to.wimp.files = example.wimp.dir,
                                  metadata = metadata)
  expect_s4_class(df.test, "phyloseq")
})

test_that("wimp_raw_to_metagenomeseq returns a MRexperiment object", {
  metadata <- read.csv(system.file(package = "epi2me2r", "extdata", "example_metadata.csv"))
  df.test <- wimp_raw_to_metagenomeseq(path.to.wimp.files = example.wimp.dir,
                                       metadata = metadata)
  expect_s4_class(df.test, "MRexperiment")
})
