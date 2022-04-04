##' Raw AMR files plus metadata to phyloseq object
##'@name amr_raw_to_metagenomeseq
##' @description given directory and metadata make phyloseq object \pkg{\link{metagenomeSeq}} package required.
##' @param path.to.amr.files path to data of raw csv files from AMRA CARD analysis
##' @param metdata data.table of metadata with "filename" and "barcode"  columns required
#' @param coveragenumber Minimum percentage of a gene that must be covered 0 to 99, default = 80
#' @param keepSNP true or false to keep AMR gene conferred by one SNP change, default = FALSE
##' @seealso \pkg{\link{metagenomeSeq}}
#' @return metagenomeSeq object for downstream analysis
#' @examples
#' \dontrun{
#' amr_raw_to_phyloseq(path.to.amr.files = path/to/amr.count.table,
#' metadata = metadata, coveragenumber = 80, keepSNP = FALSE)
#' }
#' @import data.table
#' @importFrom  phyloseq phyloseq
#' @importFrom  phyloseq otu_table
#' @importFrom phyloseq tax_table
#' @importFrom  phyloseq sample_data
#' @importFrom phyloseq phyloseq_to_metagenomeSeq
#' @importFrom Biobase AnnotatedDataFrame
#' @export

data(CARD_taxonomy, envir=environment())

amr_raw_to_phyloseq <- function(path.to.amr.files, metadata, coveragenumber=80, keepSNP=FALSE){
  #first count table
  parsed_files <- list.files(path = path.to.amr.files)
  Sample_IDs <- sub(".csv", "", parsed_files)
  i <- 1
  file_name <- paste0(parsed_files[i])
  amr.dataframe <- fread(paste0(path.to.amr.files,file_name))
  amr.dataframe <- cbind(amr.dataframe, csvname = Sample_IDs[i])

  for(i in 1:length(parsed_files)){
    file_name <- paste0(parsed_files[i])
    amr.dataframe <- fread(paste0(path.to.amr.files,file_name))
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
  amr.rawdata.reduced.subset <- amr.rawdata.reduced.subset[, list(sampleID, CVTERMID)]
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
  sampleID_names <- metadata[,"sampleID"]
  amr_count_table.t <- as.data.table(t(as.matrix(amr_count_table, rownames = "CVTERMID")), keep.rownames = "sampleID")
  dropped_mis_barcode.t <- merge(x = sampleID_names, amr_count_table.t,by = "sampleID", all.x = TRUE)
  #entire section because I could not get the df to be numeric in the count data
  dropped_mis_barcode <- as.data.frame(t(dropped_mis_barcode.t))
  dropped_col_names <- as.character(dropped_mis_barcode[1,])
  dropped_mis_barcode_table <- dropped_mis_barcode[-1, ]
  dropped_row_names <- as.character(row.names(dropped_mis_barcode_table))
  amr_table_numeric <- as.data.frame(sapply(dropped_mis_barcode_table, as.numeric))
  rownames(amr_table_numeric) <- dropped_row_names
  colnames(amr_table_numeric) <- dropped_col_names
  #taxonomy generation
  amr.CVTERMID.list <- as.data.frame(row.names(amr_table_numeric))
  setnames(amr.CVTERMID.list, "row.names(amr_table_numeric)", "CVTERMID")
  amr.CVTERMID.list$CVTERMID <- as.numeric(amr.CVTERMID.list$CVTERMID)
  CARD_taxonomy$CVTERMID <- as.numeric(CARD_taxonomy$CVTERMID)
  merged.data <- merge(x = amr.CVTERMID.list, y = CARD_taxonomy, by = "CVTERMID", all.x = TRUE)
  taxa_short <- merged.data[, c("Drug Class", "AMR Gene Family", "Resistance Mechanism", "ARO Name")]
  rownames(taxa_short) <- dropped_row_names
  #quick metadata
  metadata <- as.data.frame(metadata)
  rownames(metadata) <- metadata$sampleID
  #put it together
  OTU = otu_table(amr_table_numeric, taxa_are_rows = TRUE)
  TAX = tax_table(as.matrix(taxa_short))
  META = sample_data(metadata)
  ps <- (phyloseq(OTU, TAX, META))
  phyloseq_to_metagenomeSeq(ps)
}
