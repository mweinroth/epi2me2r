#' take raw AMR csv files to matrix
#' @param path.tofile A file path
#' @param coveragenumber A number
#' @return matrix of AMR genes at a specific coverage
#' @examples
#' read_in_amr_file("~/Desktop/my.files/", 80)
#' @import data.table

read_in_amr_files <- function(path.to.files, coveragenumber){
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
  amr.rawdata.reduced[, c('1','2', '3', '4', 'CARDontology') :=
                        do.call(Map, c(f = c, strsplit(URL, '/'))) ]
  amr.rawdata.reduced.subset <- amr.rawdata.reduced[, list(sampleID, coverage, CARDontology)]
  amr.rawdata.reduced.subset <- amr.rawdata.reduced.subset[coverage %between% c(coveragenumber, 100)]
  amr.rawdata.reduced.subset <- amr.rawdata.reduced[, list(sampleID, CARDontology)]
  mydt_wide <- suppressMessages(dcast(amr.rawdata.reduced.subset, CARDontology ~ sampleID))
  mydt_wide
}
