# epi2me2r

epi2me2r is designed to take CSV output from Oxford Nanopore's [EPI2ME](https://epi2me.nanoporetech.com/) 
and facilitate the easy import of these documents into R for downstream analysis. 
The raw data are imported either into a [phyloseq](https://bioconductor.org/packages/release/bioc/html/phyloseq.html) 
or [metagenomeSeq](https://www.bioconductor.org/packages/release/bioc/html/metagenomeSeq.html) object,
or count tables and taxonomy for other packages. 
Currently, raw data from WIMP and Antimicrobial Resistance (AMR) can be used as inputs. 

## Contents
- [Overview](#Overview)
- [Installation](#Installation)
- [Inputs](#Inputs)
- [Usage](#usage)
- [Issues](#Issues)
- [Contact](#Contact)
- [Citation](#Citation)

# Overview

There are three main types of functions in `epi2me2r`:

- **fully automated**: these functions take minimal input (just raw CSV files and metadata) and produce either `phyloseq` or `metagenomeSeq` objects for downstream analysis:
    - `raw_amr_to_phyloseq()`
        -  input: directory with all epi2me generated AMR CSV files and metadata file as a CSV
        -  output: `phyloseq` object with count table, taxonomy, and metadata
    - `raw_amr_to_metagenomeseq()`
        -  input: directory with all epi2me generated AMR CSV files and metadata file as a CSV
        -  output: `metagenomeSeq` object with count table, taxonomy, and metadata
    - `raw_wimp_to_phyloseq()`
        -  input: directory with all epi2me generated WIMP CSV files and metadata file as a CSV
        -  output: `phyloseq` object with count table, taxonomy, and metadata
    - `raw_wimp_to_metagenomeseq()`
        -  input: directory with all epi2me generated WIMP CSV files and metadata file as a CSV
        -  output: `metagenomeSeq` object with count table, taxonomy, and metadata
        
- **step-by-step**: If you are looking for just a portion of the data you can use these functions to generate only what you need:
    - `read_in_amr_files`
        -  input: directory with all epi2me generated AMR CSV files
        -  output: count matrix of all samples x AMR genes 
    - `read_in_wimp_files`
        - input: directory with all epi2me generated WIMP CSV files
        - output: count matrix of all samples x taxonomic classifications 
    - `generate_amr_taxonomy`
        - input: count table from `read_in_amr_files` or generated in the same format
        - output: ontological list of all AMR genes in the same
    - `generate_wimp_taxonomy`
        - input: count table from `read_in_wimp_files` or generated in the same format
        - output: list of of all IDed organisms in the samples by superkingdom to species

- **other**: An additional function we created but does not fit into the main workflows is `amr_read_taxonomy()`. This function reads in both AMR and WIMP raw data and adds the taxonomic information to the AMR gene if available. This function uses the paths from both the raw AMR directory containing the AMR CSVs and the WIMP directory with all the raw WIMP CSVs. The output is a taxonomic classification of all AMR genes that have both classifications available. 

## Installation

The development version of `epi2me2r` can currently be downloaded from GitHub:

```
if (!require(remotes)) install.packages('remotes')
remotes::install_github("mweinroth/epi2me2r") 
library(epi2me2r)
```

## Inputs
To use `epi2me2r` you will need your **raw data** and a **metadata file**. 

### raw data 

Raw data files are downloaded from the EPI2ME report either in the WIMP or AMR CARD tab 
(Each sample will have 2 different files if you conducted both an AMR and WIMP analysis). 
Raw data will be downloaded as a CSV file with a separate file from each run. 
If you have barcodes for your samples, multiple samples will be contained in one CSV file; 
if not, each sample will have its own file. 
**DO NOT** change the names of the files you have downloaded, 
as the file name will have the type of analysis and run number in it 
(e.g., `arma_288715.csv` [_AMR_] or `226094_1777.csv` [_WIMP_]). 

**_Place all raw data files of the same analysis type in the directory with only those files in it!_**

You will use this directory location when you import your samples. 

### metadata

The second file you will need is a metadata file describing the type of samples you have, such as sample names. 
This file may also contain other important information about your samples such as treatments.
There are four required columns if you are running both a WIMP and AMR analysis:

This file has **4** required columns that must be named as follows:

- `arma_filename` : the original amr file name without the `.csv` extension
- `arma_barcode` : the barcodes of each sample (note: if you did not barcode any of your samples, enter `none` in all of the cells). **In the AMR workflow, missing barcodes are coded as `none`**
- `wimp_filename` : the original amr file name without the `.csv` extension
- `wimp_barcode` : the barcodes of each sample (note: if you did not barcode any of your samples, enter `NA` in all of the cells). **In the WIMP workflow, missing barcodes are coded as `NA`**
- `additional information` after these four required columns, you may include any additional metadata that is important, such as treatment type, sample numbers, etc.

**A quick note about barcodes:** As talked about in the [Issues](#Issues) section below, earlier versions of the epi2me workflow list ARMA barcodes as `"barcode"` and a two digit number while WIMP barcodes are listed as `"BC"` and a two digit number. While this is no longer the case (both are written as `barcode` followed by the number), if you are using older output, see the section below on changing the barcodes to a compatible style. 

[An example CSV is available here.](https://github.com/mweinroth/epi2me2r/blob/master/inst/extdata/example_metadata.csv)

![First few rows of example CSV](https://github.com/mweinroth/epi2me2r/blob/master/inst/metadata-example.jpg)


### Usage

See [vignette](https://mweinroth.github.io/epi2me2r/articles/epi2me2r-vignette.html) for usage. 

## Issues

The current ARMA and WIMP barcodes have mostly compatible barcodes nomenclature (the only difference being ARMA files without barcodes are entered as `none` and WIMP as `NA`). However, in early epi2me versions, ARMA barcodes were listed as `"barcode"` and a two digit number while WIMP barcodes were listed as `"BC"` and a two digit number. While this is no longer the case (both are formatted as `"barcode"` followed by the number), if you are using older output, you will need to replace the `"BC"` in the WIMP files with `"barcode"`. You can do this through opening the file and `Ctrl+H` to search and replace all `"BC"` with `"barcode"` in Excel or another text editor or in R:

```
fake.data <- read.csv("226094_1777.csv") #read in data
fake.data$barcode <- sub("BC*", "barcode", fake.data$barcode) #sub
head(fake.data) #check to see if it worked
write.csv(fake.data, "226094_1777.csv") #save it (this overwrites original)
```

# Contact

Any questions or comments can be directed to Maggie Weinroth: maggie.weinroth@usda.gov

# Citation

.
