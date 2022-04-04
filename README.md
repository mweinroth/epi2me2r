# epi2me2r

epi2me2r is designed to take CSV output from Oxford Nanopore's [EPI2ME](https://epi2me.nanoporetech.com/) and facilitate the easy import of these documents into R for downstream analysis either into a [phyloseq](https://bioconductor.org/packages/release/bioc/html/phyloseq.html) or [metgenomeSeq](https://www.bioconductor.org/packages/release/bioc/html/metagenomeSeq.html) object or coutn tables and taxonomy for other packages. Currently, raw data from WIMP and AMR CARD can be used as imputs. 

## Contents
- [Overview](#Overview)
- [Installation](#Installation)
- [Inputs](#Inputs)
- [Usage](#usage)
- [Output](#output)
- [Additional Information](#Additional)
- [Contact](#Contact)

# Overview

There are three main types of functions in epi2me2r:

- **fully automated**: these functions take minimal input (usually just raw csv files and metadata) and produce either phyloseq and metagenomicSeq objects for downstream analysis:
    - `raw_amr_to_phyloseq`
    - `raw_amr_to_metagenomseq`
    - `raw_wimp_to_phyloseq`
    - `raw_wimp_to_metagenomeseq`
        
- **step-by-step**: If you are looking for just a portion of the data you can these functions to generation only what you need:
    - `read_in_amr_files`
    - `read_in_wimp_files`
    - `generate_amr_taxonomy`
    - `generate_wimp_taxonomy`

- **other**: An additional function we created but does not fit into the main workflows is `amr_read_taxonomy`. This function reads in both AMR and WIMP raw data and adds the taxonomic infromation to the AMR gene if available.  

## Installation
epi2me2r can currently be downloaded from github:
```
install.packages("devtools") 
library(devtools) 
install_github("mweinroth/epi2me2r") 
library(epi2me2r)
```

## Inputs
To use epi2me2r you will need your **raw data** and a **metadata file**. 

### raw data 
Raw data files are downloaded from the EPI2ME report:
![](https://github.com/mweinroth/epi2me2r/blob/master/screenshots.for.github/epi2me.download.report.jpg?raw=true)

Raw data will be downloaded as a csv file with a separate file from each run. If you have barcodes your samples, multiple samples will be contained in one csv file; if not each sample will have its own file. **DO NOT** change the names of the files you have downloaded as the file name will have the type of analysis and run number in it (i.e. arma_288715.csv). 

**_Place all raw data files of the same analysis type in the directory with only those files in it_**

You will use this directory location as later when you import your samples. 

### metadata
The second file you will need is a metadata file describing the type of samples you have, such as sample names. This file may also contain other important information about your samples such as treatments. There are four required columns if you are running both a WIMP and AMR analysis:

This file has **4** required columns that must been entered as seen below:

- `amr_filename` : the original amr file name without the csv extension
- `amr_barcodes` : the barcodes of each same (note if you did not barcode any of your samples enter *none* in all of the cells). **In the AMR workflow, barcodes are listed as "barcode" and a two digit number, no barcodes are entered as"none"**
- `wimp_filename` : the original amr file name without the csv extension
- `wimp_barcodes` : the barcodes of each same (note if you did not barcode any of your samples enter *NA* in all of the cells). **In the AMR workflow, barcodes are listed as "BC" and a two digit number, no barcodes are entered as"NA"**
- `additionally information` after these four required columns, you may include any additional metadata that is import, such as treatment type, sample numbers, etc.

![](https://github.com/mweinroth/epi2me2r/blob/master/screenshots.for.github/metadata.jpg)
