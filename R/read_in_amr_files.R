#' Convert raw AMR CSV files to data table
#' @name read_in_amr_files
#' @param path.to.amr.files path to folder containing raw CSV files from ARMA
#' CARD analysis
#' @param coveragenumber Minimum percentage of a gene that must be
#'  covered. Range from 0 to 99, default = 80
#' @param keepSNP TRUE or FALSE: whether to keep AMR gene conferred by one SNP
#' change, default = FALSE
#' @return data.table of AMR genes at a specific coverage with or
#' without SNP associated
#' @examples
#' \dontrun{
#' read_in_amr_files(path.to.amr.files = "~/Desktop/my.files/",
#' coveragenumber = 80, keepSNP = FALSE)
#' }
#' @import data.table
#' @export


data(CARD_taxonomy, envir=environment())

read_in_amr_files <- function(path.to.amr.files, coveragenumber=80,
                              keepSNP=FALSE) {

  # Checks for valid input. Fails with error if any are not met.
  stopifnot(coveragenumber >= 0 & coveragenumber <= 99)
  stopifnot(is.logical(keepSNP))
  stopifnot(dir.exists(path.to.amr.files))

  message(paste("Reading in raw AMR files from", path.to.amr.files))
  parsed_files <- list.files(path = path.to.amr.files)
  Sample_IDs <- sub(".csv", "", parsed_files)

  amr.rawdata <- lapply(1:length(parsed_files), function(i) {
    file_name <- paste0(parsed_files[i])
    amr.dataframe <- fread(file.path(path.to.amr.files, file_name))
    cbind(amr.dataframe, csvname = Sample_IDs[i])
  })

  amr.rawdata <- do.call(rbind, amr.rawdata)

  amr.rawdata.reduced <- amr.rawdata[, list(csvname, barcode,
                                            coverage, URL)]
  sampleidinfo <- amr.rawdata.reduced[, .(sampleID=
                                            amr.rawdata.reduced[["sampleID"]],
                                          sampleID=
                                            do.call(paste, c(.SD, sep="_"))),
                                      .SDcols= csvname:barcode]
  amr.rawdata.reduced <- cbind(sampleidinfo, amr.rawdata.reduced)
  amr.rawdata.reduced[, c('1','2', '3', '4', 'CVTERMID') :=
                        do.call(Map, c(f = c, strsplit(URL, '/'))) ]
  amr.rawdata.reduced.subset <-
    amr.rawdata.reduced[coverage %between% c(coveragenumber, 100),
                        list(sampleID, CVTERMID)]
  mydt_wide <- suppressMessages(dcast(amr.rawdata.reduced.subset,
                                      CVTERMID ~ sampleID))

  if (keepSNP) {
    mydt_wide
  } else {
    mydt_wide$CVTERMID <- as.numeric(mydt_wide$CVTERMID)
    CARD_taxonomy$CVTERMID <- as.numeric(CARD_taxonomy$CVTERMID)
    merged.data <- merge(x = mydt_wide, y = CARD_taxonomy,
                         by = "CVTERMID", all.x = TRUE)
    nosnpdata <- merged.data[mutationassociated == "no"]
    mydt_wide <- nosnpdata[, !c("ARO Accession",
                                "CARDversion", "Model Sequence ID", "Model ID",
                                "Model Name", "ARO Name",
                                "Protein Accession", "DNA Accession",
                                "AMR Gene Family", "Drug Class",
                                "Resistance Mechanism", "mutationassociated")]
    mydt_wide
  }
}
