#' Convert raw WIMP CSV files to data table
#' @param path.to.wimp.files File path to folder containing raw WIMP files
#' @return data.table of mb genes at a specific coverage
#' @examples
#' \dontrun{
#' read_in_wimp_files(path.to.wimp.files="~/Desktop/my.files")
#' }
#' @import data.table
#' @export

read_in_wimp_files <- function(path.to.wimp.files) {

  # Check inputs for validity
  stopifnot(dir.exists(path.to.wimp.files))

  message(paste("Reading in raw WIMP files from", path.to.wimp.files))
  parsed_files <- list.files(path = path.to.wimp.files)
  Sample_IDs <- sub(".csv", "", parsed_files)

  mb.rawdata <- lapply(1:length(parsed_files), function(i) {
    file_name <- paste0(parsed_files[i])
    mb.dataframe <- fread(file.path(path.to.wimp.files,file_name))
    cbind(mb.dataframe, csvname = Sample_IDs[i])
  })

  mb.rawdata <- do.call(rbind, mb.rawdata)

  total.reads <- nrow(mb.rawdata)
  mb.classified <- mb.rawdata[mb.rawdata$exit_status == "Classified",]
  classified.reads <- nrow(mb.classified)
  percentage.classified <- round((classified.reads/total.reads*100),
                                digits = 2)
  message(paste("The percentage of classified reads was",
                percentage.classified))
  mb.rawdata.reduced <- mb.rawdata[, list(csvname, barcode, taxID)]
  sampleidinfo <- mb.rawdata.reduced[, .(sampleID=
                                           mb.rawdata.reduced[["sampleID"]],
                                         sampleID=
                                           do.call(paste, c(.SD, sep="_"))),
                                     .SDcols= csvname:barcode]
  mb.rawdata.reduced <- cbind(sampleidinfo, mb.rawdata.reduced)
  mb.rawdata.reduced.subset <- mb.rawdata.reduced[, list(sampleID, taxID)]
  mydt_wide <- suppressMessages(dcast(mb.rawdata.reduced.subset,
                                        taxID ~ sampleID))
mydt_wide

}
