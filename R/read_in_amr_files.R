#' take raw AMR csv files to matrix
#' @name read_in_amr_files
#' @param path.to.file A file path to AMR raw data csv
#' @param coveragenumber Minimum percentage of a gene that must be covered 0 to 99
#' @param keepSNP true or false to keep AMR gene conferred by one SNP change
#' @return matrix of AMR genes at a specific coverage with or without SNP associated
#' @examples
#' \dontrun{
#' read_in_amr_files(path.to.file = "~/Desktop/my.files/", coveragenumber = 80, keepSNP = FALSE)
#' }
#' @import data.table
#' @export


data(CARD_taxonomy, envir=environment())

read_in_amr_files <- function(path.to.files, coveragenumber=80, keepSNP=FALSE){
  parsed_files <- list.files(path = path.to.files)
  Sample_IDs <- sub(".csv", "", parsed_files)
  i <- 1
  file_name <- paste0(parsed_files[i])
  amr.dataframe <- fread(paste0(path.to.files,file_name))
  amr.dataframe <- cbind(amr.dataframe, csvname = Sample_IDs[i])

  for(i in 1:length(parsed_files)){
    file_name <- paste0(parsed_files[i])
    amr.dataframe <- fread(paste0(path.to.files,file_name))
    amr.dataframe <- cbind(amr.dataframe, csvname = Sample_IDs[i])
    if(i == 1){
      amr.rawdata <- amr.dataframe
    } else {
      amr.rawdata <- rbind(amr.rawdata, amr.dataframe)
    }
  }
  amr.rawdata.reduced <- amr.rawdata[, list(csvname, barcode, coverage, URL)]
  sampleidinfo <- amr.rawdata.reduced[, .(sampleID=amr.rawdata.reduced[["sampleID"]],
                                          sampleID=do.call(paste, c(.SD, sep="_"))),
                                      .SDcols= csvname:barcode]
  amr.rawdata.reduced <- cbind(sampleidinfo, amr.rawdata.reduced)
  amr.rawdata.reduced[, c('1','2', '3', '4', 'CVTERMID') :=
                        do.call(Map, c(f = c, strsplit(URL, '/'))) ]
  amr.rawdata.reduced.subset <- amr.rawdata.reduced[, list(sampleID, coverage, CVTERMID)]
  amr.rawdata.reduced.subset <- amr.rawdata.reduced.subset[coverage %between% c(coveragenumber, 100)]
  amr.rawdata.reduced.subset <- amr.rawdata.reduced[, list(sampleID, CVTERMID)]
  mydt_wide <- suppressMessages(dcast(amr.rawdata.reduced.subset, CVTERMID ~ sampleID))
  {
    if (keepSNP == TRUE) {
      print(mydt_wide)
    } else if (keepSNP == FALSE) {
      mydt_wide$CVTERMID <- as.numeric(mydt_wide$CVTERMID)
      CARD_taxonomy$CVTERMID <- as.numeric(CARD_taxonomy$CVTERMID)
      merged.data <- merge(x = mydt_wide, y = CARD_taxonomy, by = "CVTERMID", all.x = TRUE)
      nosnpdata <- merged.data[mutationassociated == "no"]
      count_nosnpdata <- nosnpdata[, !c("ARO Accession", "CARDversion", "Model Sequence ID", "Model ID",
                                        "Model Name", "ARO Name", "Protein Accession", "DNA Accession",
                                        "AMR Gene Family", "Drug Class", "Resistance Mechanism", "mutationassociated")]
      print(count_nosnpdata)
    }
  }
}
