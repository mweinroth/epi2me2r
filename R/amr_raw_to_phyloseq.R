##' Raw AMR files plus metadata to phyloseq object
##'@name amr_raw_to_phyloseq
##' @description given directory and metadata make phyloseq object \pkg{\link{phyloseq}} package required.
##' @param x Data.
##' @param y More data.
##' @seealso \pkg{\link{phyloseq}}
##' @export
#' @param amr.count.table count of AMR genes with CVTERMID in the first row and samples on columns
#' @param verbose true or false to keep all columns of AMR gene information or not
#' @return data table of taxonomy of AMR genes in count table
#' @examples
#' \dontrun{
#' generate_amr_taxonomy(amr.count.table = amr.count.table, verbose = FALSE)
#' }
#' @import data.table
#' @export

data(CARD_taxonomy, envir=environment())



data(CARD_taxonomy, envir=environment())

amr_raw_to_phyloseq <- function(path.to.files, metadata, coveragenumber = 80, keepSNP = FALSE){
  #first count table
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
      amr_count_table <- print(mydt_wide)
    } else if (keepSNP == FALSE) {
      mydt_wide$CVTERMID <- as.numeric(mydt_wide$CVTERMID)
      CARD_taxonomy$CVTERMID <- as.numeric(CARD_taxonomy$CVTERMID)
      merged.data <- merge(x = mydt_wide, y = CARD_taxonomy, by = "CVTERMID", all.x = TRUE)
      nosnpdata <- merged.data[mutationassociated == "no"]
      count_nosnpdata <- nosnpdata[, !c("ARO Accession", "CARDversion", "Model Sequence ID", "Model ID",
                                        "Model Name", "ARO Name", "Protein Accession", "DNA Accession",
                                        "AMR Gene Family", "Drug Class", "Resistance Mechanism", "mutationassociated")]
      amr_count_table <- count_nosnpdata
    }
  }
  #remove mis-barcoded samples
  metadata$sampleID <- paste(metadata$filename, metadata$barcode, sep = "_")
  # sampleID_names <- metadata[,"sampleID"]
  # amr_count_table.t <- as.data.table(t(as.matrix(csvtry2, rownames = "CVTERMID")), keep.rownames = "sampleID")
  #dropped_mis_barcode.t <- merge(x = sampleID_names, amr_count_table.t,by = "sampleID", all.x = TRUE)
  # dropped_mis_barcode <- as.data.table(t(dropped_mis_barcode.t), keep.rownames = TRUE)
  #colnames(dropped_mis_barcode) <- dropped_mis_barcode[1,]
  # dropped_mis_barcode_table <- dropped_mis_barcode[-1, ]

  #taxonomy generation
  amr.CVTERMID.list <- amr_count_table[,1]
  amr.CVTERMID.list$CVTERMID <- as.numeric(amr.CVTERMID.list$CVTERMID)
  CARD_taxonomy$CVTERMID <- as.numeric(CARD_taxonomy$CVTERMID)
  merged.data <- merge(x = amr.CVTERMID.list, y = CARD_taxonomy, by = "CVTERMID", all.x = TRUE)
  taxa_short <- merged.data[, c("CVTERMID","Drug Class", "AMR Gene Family", "Resistance Mechanism", "ARO Name")]

  #put it together
  OTU = otu_table(amr_count_table, taxa_are_rows = TRUE)
  TAX = tax_table(as.matrix(taxa_short))
  META = sample_data(metadata)
  phyloseq(OTU, TAX)

}

#phyloseq_object <- amr_raw_to_phyloseq(path.to.files = "raw-amr/", metadata = metadata)

#metadata <- read_csv("meta-amr.csv")
