# epi2me2r

epi2me2r is designed to take CSV output from Oxford Nanopore's [EPI2ME](https://epi2me.nanoporetech.com/) and facilitate the easy import of these documents into R for downstream analysis. Currently, raw data from WIMP and AMR CARD can be used. 

## Contents
- [Installation](#Installation)
- [Inputs](#Inputs)
- [Installation](#installation)
- [Usage](#usage)
- [Output](#output)


## Installation
epi2me2r can currently be downloaded from github:
```
install.packages("devtools") 
library(devtools) 
install_github("mweinroth/epi2me2") 
library(epi2me2)
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
The second file you will need is a metadata file describing the type of samples you have, such as sample names. This file may also contain other important information about your samples such as treatments.
