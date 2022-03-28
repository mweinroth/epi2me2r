#' take raw mb csv files to matrix
#' @param path.tofile A file path
#' @return matrix of mb genes at a specific coverage
#' @examples
#' \dontrun{
#' read_in_mb_file(path.to.files="~/Desktop/my.files")
#' }
#' @import data.table
#' @export

read_in_wimp_files <- function(path.to.files){
  parsed_files <- list.files(path = path.to.files)
  Sample_IDs <- sub(".csv", "", parsed_files)
  i <- 1
  file_name <- paste0(parsed_files[i])
  mb.dataframe <- fread(paste0(path.to.files,file_name))
  mb.dataframe <- cbind(mb.dataframe, csvname = Sample_IDs[i])
  for(i in 1:length(parsed_files)){
    file_name <- paste0(parsed_files[i])
    mb.dataframe <- fread(paste0(path.to.files,file_name))
    mb.dataframe <- cbind(mb.dataframe, csvname = Sample_IDs[i])
    if(i == 1){
      mb.rawdata <- mb.dataframe
    } else {
      mb.rawdata <- rbind(mb.rawdata, mb.dataframe)
    }
  }

  total.reads <- nrow(mb.rawdata)
  mb.classified <- mb.rawdata[mb.rawdata$exit_status == "Classified",]
  classified.reads <- nrow(mb.classified)
  percentage.classifed <- round((classified.reads/total.reads*100), digits = 2)
  message(paste("The percentage of classifed reads was", percentage.classifed))
  mb.rawdata.reduced <- mb.rawdata[, list(csvname, barcode, taxID)]
  sampleidinfo <- mb.rawdata.reduced[, .(sampleID=mb.rawdata.reduced[["sampleID"]],
                                         sampleID=do.call(paste, c(.SD, sep="_"))),
                                     .SDcols= csvname:barcode]
  mb.rawdata.reduced <- cbind(sampleidinfo, mb.rawdata.reduced)
  mb.rawdata.reduced.subset <- mb.rawdata.reduced[, list(sampleID, taxID)]
  mydt_wide <- suppressMessages(dcast(mb.rawdata.reduced.subset, taxID ~ sampleID))
  mydt_wide
}
