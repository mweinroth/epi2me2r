##' Raw AMR files plus metadata to phyloseq object
##'@name amr_raw_to_metagenomeseq
##' @description Given directory and metadata, make phyloseq object.
##' \pkg{\link{metagenomeSeq}} package required.
##' @param path.to.amr.files path to data of raw csv files from AMRA
##' CARD analysis
##' @param metadata data.table of metadata with "filename" and "barcode"
##'  columns required
#' @param coveragenumber Minimum percentage of a gene that must be
#'  covered. Range from 0 to 99, default = 80
#' @param keepSNP TRUE or FALSE: whether to keep AMR gene conferred by one SNP
#' change, default = FALSE
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
#' @import metagenomeSeq
#' @export

data(CARD_taxonomy, envir=environment())

amr_raw_to_metagenomeseq <- function(path.to.amr.files, metadata,
                                     coveragenumber=80, keepSNP=FALSE){

  # Checks for valid input. Fails with error if any are not met.
  stopifnot(coveragenumber >= 0 & coveragenumber <= 99)
  stopifnot(is.logical(keepSNP))
  stopifnot(dir.exists(path.to.amr.files))
  stopifnot(is.data.frame(metadata))

  if (any(!c('filename', 'barcode') %in% names(metadata))) {
    stop('metadata does not have columns named "filename" and "barcode".')
  }

  amr_count_table <- read_in_amr_files(path.to.amr.files, coveragenumber, keepSNP)

  #remove mis-barcoded samples
  metadata$sampleID <- paste(metadata$arma_filename,
                             metadata$amra_barcode, sep = "_")
  sampleID_names <- as.data.frame(metadata$sampleID)
  colnames(sampleID_names) <- "sampleID"
  amr_count_table.t <- as.data.table(t(as.matrix(amr_count_table,
                                                 rownames = "CVTERMID")),
                                     keep.rownames = "sampleID")
  dropped_mis_barcode.t <- merge(x = sampleID_names, amr_count_table.t,
                                 by = "sampleID", all.x = TRUE)
  #entire section because I could not get the df to be numeric in the count data
  dropped_mis_barcode <- as.data.frame(t(dropped_mis_barcode.t))
  dropped_col_names <- as.character(dropped_mis_barcode[1,])
  dropped_mis_barcode_table <- dropped_mis_barcode[-1, ]
  dropped_row_names <- as.character(row.names(dropped_mis_barcode_table))
  amr_table_numeric <- as.data.frame(sapply(dropped_mis_barcode_table,
                                            as.numeric))
  rownames(amr_table_numeric) <- dropped_row_names
  colnames(amr_table_numeric) <- dropped_col_names
  #taxonomy generation
  amr.CVTERMID.list <- as.data.frame(row.names(amr_table_numeric))
  setnames(amr.CVTERMID.list, "row.names(amr_table_numeric)", "CVTERMID")
  amr.CVTERMID.list$CVTERMID <- as.numeric(amr.CVTERMID.list$CVTERMID)
  CARD_taxonomy$CVTERMID <- as.numeric(CARD_taxonomy$CVTERMID)
  merged.data <- merge(x = amr.CVTERMID.list,
                       y = CARD_taxonomy, by = "CVTERMID", all.x = TRUE)
  taxa_short <- merged.data[, c("Drug Class",
                                "AMR Gene Family", "Resistance Mechanism",
                                "ARO Name")]
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
