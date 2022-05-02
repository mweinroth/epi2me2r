library(epi2me2r)
test_that("amr_raw_to_phyloseq gives phyloseq object", {
  metadata <- read.csv("testdata/example_metadata.csv")
  df.test <- amr_raw_to_phyloseq(path.to.amr.files = "testdata/example_amr_data/",
                                 metadata = metadata)
  expect_s4_class(df.test, "phyloseq")
})
