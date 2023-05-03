#' assign taxonomy for phylogenetic and AMR for each read
#'@name amr_read_taxonomy
#' @description Given raw data for AMR and WIMP, provides full AMR and
#' taxon info for those reads that assign to both
#' @param path.to.wimp.files path to folder containing raw CSV files from WIMP analysis
#' @param path.to.amr.files path to folder containing raw CSV files from ARMA analysis
#' @param coveragenumber Minimum percentage of a gene that must be
#'  covered. Range from 0 to 99, default = 80
#' @return data.frame with double classified reads
#' @examples
#' \dontrun{
#' amr_read_taxonomy(path.to.wimp.files = "path/to/wimp_files",
#' path.to.amr.files = "path/to/amr_files", coveragenumber = 80)
#' }
#' @import data.table
#' @import stats
#' @export

data(CARD_taxonomy, envir=environment())

amr_read_taxonomy <- function(path.to.wimp.files, path.to.amr.files,
                              coveragenumber=80){

  # Checks for valid input. Fails if any is invalid.
  stopifnot(coveragenumber >= 0 & coveragenumber <= 99)
  stopifnot(dir.exists(path.to.wimp.files))
  stopifnot(dir.exists(path.to.amr.files))

  read_in_amr_readid <- function(path.to.amr.directory, coveragepercentage) {
    parsed_files <- list.files(path = path.to.amr.directory)
    Sample_IDs <- sub(".csv", "", parsed_files)
    i <- 1
    file_name <- paste0(parsed_files[i])
    amr.dataframe <- fread(file.path(path.to.amr.directory,file_name))
    amr.dataframe <- cbind(amr.dataframe, csvname = Sample_IDs[i])

    for(i in 1:length(parsed_files)){
      file_name <- paste0(parsed_files[i])
      amr.dataframe <- fread(file.path(path.to.amr.directory,file_name))
      colnames(amr.dataframe) <- tolower(colnames(amr.dataframe))
      amr.dataframe <- cbind(amr.dataframe, csvname = Sample_IDs[i])
      if(i == 1){
        amr.rawdata <- amr.dataframe
      } else {
        amr.rawdata <- rbind(amr.rawdata, amr.dataframe)
      }
    }
    amr.rawdata.reduced <- amr.rawdata[, list(csvname, barcode,
                                              coverage, url, read_id)]
    sampleidinfo <- amr.rawdata.reduced[, .(sampleid=
                                              amr.rawdata.reduced[["sampleid"]],
                                            sampleid=
                                              do.call(paste, c(.SD, sep="_"))),
                                        .SDcols= csvname:barcode]
    amr.rawdata.reduced <- cbind(sampleidinfo, amr.rawdata.reduced)
    amr.rawdata.reduced[, c('1','2', '3', '4', 'CVTERMID') :=
                          do.call(Map, c(f = c, strsplit(url, '/'))) ]
    amr.rawdata.reduced.subset <- amr.rawdata.reduced[, list(sampleid,
                                                             coverage,
                                                             CVTERMID,
                                                             read_id)]
    amr.rawdata.reduced.subset <- amr.rawdata.reduced.subset[
      coverage %between% c(coveragepercentage, 100)]
    amr.rawdata.reduced.subset <- amr.rawdata.reduced.subset[,
                                                             list(read_id,
                                                                  sampleid,
                                                                  CVTERMID)]
  }
  read_in_wimp_files <- function(path.to.wimp.directory) {
    message(paste("Reading in raw files from", path.to.wimp.directory))
    parsed_files <- list.files(path = path.to.wimp.directory)
    Sample_IDs <- sub(".csv", "", parsed_files)

    mb.rawdata <- lapply(1:length(parsed_files), function(i) {
      file_name <- paste0(parsed_files[i])
      mb.dataframe <- fread(file.path(path.to.wimp.files,file_name))
      # Change all col names to lower, avoiding mismatches
      colnames(mb.dataframe) <- tolower(colnames(mb.dataframe))
      # get only necessary information for wimp
      # new reports generate other columns like species_ID and genus_ID
      mb.dataframe =  dplyr:::select(mb.dataframe,filename,readid,runid,exit_status,barcode,taxid)
      cbind(mb.dataframe, csvname = Sample_IDs[i])
    })

    mb.rawdata <- do.call(rbind, mb.rawdata)

    total.reads <- nrow(mb.rawdata)
    # Some outputs include 'Classification successful'
    # again to lower to avoid mismatches
    mb.classified <- mb.rawdata[tolower(mb.rawdata$exit_status) %in% c("classified", "classification successful"),]
    classified.reads <- nrow(mb.classified)
    mb.rawdata.reduced <- mb.rawdata[, list(csvname, barcode, taxid, readid)]
    sampleidinfo <- mb.rawdata.reduced[, .(sampleid=
                                             mb.rawdata.reduced[["sampleid"]],
                                           sampleid=
                                             do.call(paste, c(.SD, sep="_"))),
                                       .SDcols= csvname:readid]
    mb.rawdata.reduced <- cbind(sampleidinfo, mb.rawdata.reduced)
    mb.rawdata.reduced.subset <- mb.rawdata.reduced[, list(readid, taxid)]
    setnames(mb.rawdata.reduced.subset, "readid", "read_id")
  }
  #read in
  amr.file <- read_in_amr_readid(path.to.amr.directory = path.to.amr.files,
                                 coveragepercentage=coveragenumber)
  wimp.file <- read_in_wimp_files(path.to.wimp.directory = path.to.wimp.files)

  combo_wimp_amr <- merge(x = amr.file,
                          y = wimp.file,
                          by = "read_id", all.x = TRUE)
  combo_classified_only <- na.omit(combo_wimp_amr)

  combo_classified_only$CVTERMID <- as.numeric(combo_classified_only$CVTERMID)
  CARD_taxonomy$CVTERMID <- as.numeric(CARD_taxonomy$CVTERMID)
  merged.data <- merge(x = combo_classified_only, y = CARD_taxonomy,
                       by = "CVTERMID", all.x = TRUE)
  taxa_short <- merged.data[, c("read_id", "CVTERMID",
                                "Drug Class", "AMR Gene Family",
                                "Resistance Mechanism", "ARO Name",
                                "mutationassociated", "taxid")]

  #WIMP taxonomy
  mb.taxonIDneeded <- as.numeric(taxa_short$taxid)
  message("Now downloading and putting together the NCBI database. This might take a while...")
  prepareDatabase(getAccessions=FALSE, indexTaxa=TRUE)
  message("Assigning all taxID in count matrix to fill taxonomy. You might want to take a break...")
  full.taxon.wimp <- getTaxonomy(mb.taxonIDneeded,'nameNode.sqlite')
  full.taxon.wimp.dt <- as.data.table(full.taxon.wimp,
                                      keep.rownames = "taxid")
  full.taxon.wimp.dt$taxid <- as.numeric(full.taxon.wimp.dt$taxid)
  merged.data.both <- merge(x = taxa_short,
                            y = full.taxon.wimp.dt,
                            by = "taxid", all = TRUE, allow.cartesian=TRUE)
  merged.data.both <- unique.data.frame(merged.data.both)
  merged.data.both <- merged.data.both[, c("read_id", "CVTERMID","Drug Class",
                                           "AMR Gene Family",
                                           "Resistance Mechanism",
                                           "ARO Name", "mutationassociated",
                                           "taxid", "superkingdom",
                                           "phylum", "class", "order", "family",
                                           "genus", "species")]
  merged.data.both
}
