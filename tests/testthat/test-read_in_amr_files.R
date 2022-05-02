library(epi2me2r)
system.file(package="epi2me2r")
test_that("read_in_amr works", {
  df.test <- read_in_amr_files(path.to.amr.files = "testdata/example_amr_data/")
  expect_s3_class(df.test, "data.frame")
})
