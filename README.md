# epi2me2r

epi2me2r is designed to take CSV output from Oxford Nanopore's [EPI2ME](https://epi2me.nanoporetech.com/) and facilitate the easy import of these documents into R for downstream analysis either into a [phyloseq](https://bioconductor.org/packages/release/bioc/html/phyloseq.html) or [metgenomeSeq](https://www.bioconductor.org/packages/release/bioc/html/metagenomeSeq.html) object or count tables and taxonomy for other packages. Currently, raw data from WIMP and Antimicrobial Resistance (AMR) can be used as inputs. 

## Contents
- [Overview](#Overview)
- [Installation](#Installation)
- [Inputs](#Inputs)
- [Usage](#usage)
- [Options](#Options)
- [Issues](#Issues)
- [Contact](#Contact)
- [Citation](#Citation)

# Overview

There are three main types of functions in epi2me2r:

- **fully automated**: these functions take minimal input (just raw csv files and metadata) and produce either phyloseq and metagenomicSeq objects for downstream analysis:
    - `raw_amr_to_phyloseq`
        -  input: directory with all epi2me generated AMR csv files and metadata file as a csv
        -  output: phyloseq object with count table, taxonomy, and metadata
    - `raw_amr_to_metagenomseq`
        -  input: directory with all epi2me generated AMR csv files and metadata file as a csv
        -  output: metagenomeseq object with count table, taxonomy, and metadata
    - `raw_wimp_to_phyloseq`
        -  input: directory with all epi2me generated WIMP csv files and metadata file as a csv
        -  output: phyloseq object with count table, taxonomy, and metadata
    - `raw_wimp_to_metagenomeseq`
        -  input: directory with all epi2me generated WIMP csv files and metadata file as a csv
        -  output: metagenomeseq object with count table, taxonomy, and metadata
        
- **step-by-step**: If you are looking for just a portion of the data you can these functions to generation only what you need:
    - `read_in_amr_files`
        -  input: directory with all epi2me generated AMR csv files
        -  output: count matrix of all samples x AMR genes 
    - `read_in_wimp_files`
        - input: directory with all epi2me generated WIMP csv files
        - output: count matrix of all samples x taxinomic classifications 
    - `generate_amr_taxonomy`
        - input: count table from `read_in_amr_files` or generated in the same format
        - output: ontological list of all AMR genes in the same
    - `generate_wimp_taxonomy`
        - input: count table from `read_in_wimp_files` or generated in the same format
        - output: list of of all IDed organisms in the samples by superkingdom to species

- **other**: An additional function we created but does not fit into the main workflows is `amr_read_taxonomy`. This function reads in both AMR and WIMP raw data and adds the taxonomic information to the AMR gene if available. This function uses the paths from both the raw AMR directory containing the AMR  csv and the WIMP directory with all the raw WIMP csv. The Output is a taxonomic classification of all AMR genes that have both classifications available. 

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
Raw data files are downloaded from the EPI2ME report either in the WIMP or AMR CARD tab (Each sample will have 2 different files if you conducted both an AMR and WIMP analysis). Raw data will be downloaded as a csv file with a separate file from each run. If you have barcodes your samples, multiple samples will be contained in one csv file; if not each sample will have its own file. **DO NOT** change the names of the files you have downloaded as the file name will have the type of analysis and run number in it (i.e. arma_288715.csv [_AMR_] or 226094_1777.csv [_WIMP_]). 

**_Place all raw data files of the same analysis type in the directory with only those files in it_**

You will use this directory location when you import your samples. 

### metadata
The second file you will need is a metadata file describing the type of samples you have, such as sample names. This file may also contain other important information about your samples such as treatments. There are four required columns if you are running both a WIMP and AMR analysis:

This file has **4** required columns that must been entered as seen below:

- `amr_filename` : the original amr file name without the csv extension
- `amr_barcodes` : the barcodes of each same (note if you did not barcode any of your samples enter *none* in all of the cells). **In the AMR workflow, barcodes are listed as "barcode" and a two digit number, no barcodes are entered as "none"**
- `wimp_filename` : the original amr file name without the csv extension
- `wimp_barcodes` : the barcodes of each same (note if you did not barcode any of your samples enter *NA* in all of the cells). **In the AMR workflow, barcodes are listed as "BC" and a two digit number, no barcodes are entered as "NA"**
- `additionally information` after these four required columns, you may include any additional metadata that is import, such as treatment type, sample numbers, etc.

An example csv is avalible [here](https://github.com/mweinroth/epi2me2r/blob/master/data/example_metadata.csv)

![](https://github.com/mweinroth/epi2me2r/blob/master/inst/metadata-example.jpg)


### Usage


### Full Automated data import

For both AMR and WIMP data, the raw csv's downloaded from the epi2mw website need to be in their own directory (without any other files). Note that if you are processing both WIMP and ARM data you will need **two** directories one for each set of data. 

#### AMR data

##### amr_raw_to_phyloseq

The read in for AMR data requires a directory and a metadata file. The directory will be one that only has the csv files generated from EPI2ME in it. And example of the metadata file is above. The data we will be using is from an example run on the EPI2ME pipeline. There are four options:

- `path.to.amr.files` *required*: the path to the raw csv files (i.e. "Desktop/raw_data/")
- `metadata` *required*: metadata following format above as a data.frame
- `coveragenumber` *optional* : the total length of the gene that must be present for it to be included in the count table, the default is 80%; this argument takes any number from 1 to 99. Default is `80`.
- `keepSNP` *optional* : if gene that are considered resistance only with a SNP mutation. Default is `FALSE` (does not include these genes) but can be changed to `TRUE` to include these genes. 


```{r amr phylo}
ps.amr.object <- amr_raw_to_phyloseq(path.to.amr.files = "../data/example_amr_data/",
                                     metadata = epi2me.metadata,
                                     coveragenumber = 80, 
                                     keepSNP = FALSE)
```

##### amr_raw_to_metagenomeseq

The same parameters for importing to metagenomseq are used as those that were used above in the `amr_raw_to_phyloseq` function:

```{r amr mgs}
mgs.amr.object <- amr_raw_to_metagenomeseq(path.to.amr.files = "../data/example_amr_data/",
                                     metadata = epi2me.metadata,
                                     coveragenumber = 80, 
                                     keepSNP = FALSE)
mgs.amr.object
```


#### WIMP

##### wimp_raw_to_phyloseq

WIMP files used the same idea as the AMR files but use the package [taxonomizr](https://github.com/sherrillmix/taxonomizr) to help add taxonomic hierarchical information. 

The read in for WIMP data requires a directory and a metadata file. The directory will be one that only has the csv files generated from EPI2ME in it. And example of the metadata file is above. The data we will be using is from an example run on the EPI2ME pipeline. There are four options:

- `path.to.wimp.files` *required*: the path to the raw csv files (i.e. "Desktop/raw_data/")
- `metadata` *required*: metadata following format above as a data.frame
- `keep.unclassifed` *optional* : keep genes that do not classify or do not classify beyond a superkingdom. Default is `FALSE` (does not include these reads) but can be changed to `TRUE` to include these reads 
- `keep.human` *optional* : Keep reads associated with homo sapien (for microbiome usually considered a contaminant) Default is `FALSE` (does not include human associated reads) but can be changed to `TRUE` to include these reads. 

```{r wimp phylo, eval = FALSE}
ps.wimp.object <- wimp_raw_to_phyloseq(path.to.wimp.files = "../data/example_wimp_data/",
                                     metadata = epi2me.metadata,
                                     keep.unclassifed = FALSE, 
                                     keep.human = FALSE)
```                                     

##### wimp_raw_to_metagenomeseq

Like the functions for AMR, the same parameters for importing to metagenomseq are used as those that were used above in the `wimp_raw_to_phyloseq` function:

```{r wimp mgs, eval = FALSE}
mgs.wimp.object <- wimp_raw_to_metagenomeseq(path.to.wimp.files = "../data/example_wimp_data/",
                                     metadata = epi2me.metadata,
                                     keep.unclassifed = FALSE, 
                                     keep.human = FALSE)
```

### Step-by-step import

In some cases you might not want a phyloseq or metagenomeseq object but instead just a count matrix or taxonomic list. In these cases you can use the below functions. 

#### AMR data

##### read_in_amr_file

This takes the directory the AMR csv files are in and creates a count matrix that can be used in down stream analysis. The inputs are similar to those in the previous examples (but does not require metadata):

- `path.to.amr.files` *required*: the path to the raw csv files (i.e. "Desktop/raw_data/")
- `coveragenumber` *optional* : the total length of the gene that must be present for it to be included in the count table, the default is 80%; this argument takes any number from 1 to 99. Default is `80`.
- `keepSNP` *optional* : if gene that are considered resistance only with a SNP mutation. Default is `FALSE` (does not include these genes) but can be changed to `TRUE` to include these genes. 

```{r amr read in}
amr.count.table <- read_in_amr_files(path.to.amr.files = "../data/example_amr_data/",
                                     coveragenumber = 80, 
                                     keepSNP = FALSE)
```

##### generate_amr_taxonomy

This function assigns AMR taxonomic hierarchical information from [CARD](https://card.mcmaster.ca/) using a count table with CV TERM ID's as the first column (CVTERMID). Only one input is needed:

- `amr.count.table` *required*: data frame of generated with `amr.count.table` or that has  `CVTERMID` as the first column for AMR taxonomic assignment
- `verbose` *optional* : only a subset of coumn names are included in the output by default: (`CVTERMID`,`Drug Class`, `AMR Gene Family`, `Resistance Mechanism`, and `ARO Name`) if `verbose=TRUE` 13 columns are returned. 

```{r amr tax read in}
amr.taxonomy <- generate_amr_taxonomy(amr.count.table = amr.count.table,
                                         verbose = FALSE)
```

#### WIMP data

##### read_in_wimp_file

This takes the directory the WIMP csv files are in and creates a count matrix that can be used in down stream analysis. The inputs are similar to those in the previous examples (but does not require metadata):

- `path.to.wimp.files` *required*: the path to the raw csv files (i.e. "Desktop/raw_data/")

```{r wimp read in}
wimp.count.table <- read_in_wimp_files(path.to.wimp.files = "../data/example_wimp_data/")
head(wimp.count.table)
```

##### generate_wimp_taxonomy

This function assigns phylogenetic taxonomic hierarchical information with the help of `taxonomizr` a count table with NCBI taxonomic ID's (`taxID` as the first column is required.

- `wimp.count.table` *required*: data frame of generated with `wimp.count.table` or that has  `taxID` as the first column for phylogenetic taxonomic assignment

```{r wimp tax, eval=F}
wimp.taxonomy <- generate_wimp_taxonomy(wimp.count.table = wimp.count.table)
```

## Options 

Some of the fucntions have additional options that are set to a default setting put can be changed. 

## Issues

metadata barcodes between versions of epi2me

# Contact

Any questions or comments can be directed to Maggie Weinroth: maggie.weinroth@usda.gov

# Citation

.
