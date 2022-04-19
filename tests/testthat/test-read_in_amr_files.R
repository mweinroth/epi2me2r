test_that("read_in_amr works", {
  df.test <- read_in_amr_files(path.to.amr.files = "../../inst/example_amr_data/")
  expect_s3_class(df.test, "data.frame")
})
