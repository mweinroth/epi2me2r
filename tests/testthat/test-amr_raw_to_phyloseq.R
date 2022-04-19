test_that("amr_raw_to_phyloseq gives phyloseq object", {
  metadata <- read.csv("../../inst/example_metadata.csv")
  df.test <- amr_raw_to_phyloseq(path.to.amr.files = "../../inst/example_amr_data/",
                                 metadata = metadata)
  expect_s4_class(df.test, "phyloseq")
})
