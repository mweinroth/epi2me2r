# unload previous downloaded epi2me2r
unloadNamespace("epi2me2r")

# reload required libraries
library(data.table)
library(testthat)
library(Biobase)
library(taxonomizr)
library(tidyverse)
library(phyloseq)


# Go to source file
path <- rstudioapi::getActiveDocumentContext()$path
Encoding(path) <- "UTF-8"
setwd(dirname(path))
# return to home pwd
setwd("../../../")

# Manually load modified scripts
source("epi2me2r/R/globals.R")
source("epi2me2r/R/data.R")
source("epi2me2r/R/read_in_wimp_files.R")
source("epi2me2r/R/wimp_raw_to_metagenomeseq.R")
source("epi2me2r/R/wimp_raw_to_phyloseq.R")
source("epi2me2r/R/generate_wimp_taxonomy.R")

# define variables
example.amr.dir <- "epi2me2r/inst/extdata/example_amr_data/"
example.wimp.dir <- "epi2me2r/inst/extdata/example_wimp_data/"

# load workspace CARD
load("epi2me2r/data/CARD_taxonomy.RData")

# Perform tests
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
  metadata <- read.csv("epi2me2r/inst/extdata/example_metadata.csv")
  df.test <- wimp_raw_to_phyloseq(path.to.wimp.files = example.wimp.dir,
                                  metadata = metadata)
  expect_s4_class(df.test, "phyloseq")
})

test_that("wimp_raw_to_metagenomeseq returns a MRexperiment object", {
  metadata <- read.csv("epi2me2r/inst/extdata/example_metadata.csv")
  df.test <- wimp_raw_to_metagenomeseq(path.to.wimp.files = example.wimp.dir,
                                       metadata = metadata)
  expect_s4_class(df.test, "MRexperiment")
})
