##' Raw WIMP files plus metadata to phyloseq object
##'@name wimp_raw_to_phyloseq
##' @description Given wimp directory and metadata, make phyloseq object.
##' \pkg{\link{phyloseq}} package required.
##' @param path.to.wimp.files path to data of raw csv files from WIMP
##' analysis
##' @param metadata dataframe of metadata with "filename" and "barcode"
##'  columns required
#' @param keep.unclassified TRUE or FALSE: whether to keep reads that do not
#' classify below phylum, default = FALSE
#' @param keep.human TRUE or FALSE: whether to keep reads that are classified as
#' human, default = FALSE
##' @seealso \pkg{\link{phyloseq}}
#' @return phyloseq object for downstream analysis with WIMP data
#' @examples
#' \dontrun{
#' wimp_raw_to_phyloseq(path.to.wimp.files = path/to/wimpfiles,
#' metadata = metadata, keep.unclassified = FALSE, keep.human = FALSE)
#' }
#' @import data.table
#' @import taxonomizr
#' @importFrom  phyloseq phyloseq
#' @importFrom  phyloseq otu_table
#' @importFrom phyloseq tax_table
#' @importFrom  phyloseq sample_data
#' @export


wimp_raw_to_phyloseq <- function(path.to.wimp.files, metadata,
                                 keep.unclassified=FALSE, keep.human=FALSE){

  # Checks for valid input. Fails with error if any are not met.
  stopifnot(dir.exists(path.to.wimp.files))
  stopifnot(is.data.frame(metadata))
  stopifnot(is.logical(keep.unclassified))
  stopifnot(is.logical(keep.human))

  if (any(!c('filename', 'barcode') %in% names(metadata))) {
    stop('metadata does not have columns named "filename" and "barcode".')
  }

  wimp_count_table <- read_in_wimp_files(path.to.wimp.files)

  #remove barcodes not in metadata
  mb.metadata <- metadata
  mb.metadata$sampleID <- paste(mb.metadata$wimp_filename,
                                mb.metadata$wimp_barcode, sep = "_")
  sampleID_names <- as.data.frame(mb.metadata$sampleID)
  colnames(sampleID_names) <- "sampleID"
  mb_count_table.t <- as.data.table(t(as.matrix(wimp_count_table,
                                                rownames = "taxID")),
                                    keep.rownames = "sampleID")
  dropped_mis_barcode.t <- merge(x = sampleID_names, mb_count_table.t,
                                 by = "sampleID", all.x = TRUE)
  #entire section because I could not get the df to be numeric in count data
  dropped_mis_barcode <- as.data.frame(t(dropped_mis_barcode.t))
  dropped_col_names <- as.character(dropped_mis_barcode[1,])
  dropped_mis_barcode_table <- dropped_mis_barcode[-1, ]
  dropped_row_names <- as.character(row.names(dropped_mis_barcode_table))
  mb_table_numeric <- as.data.frame(sapply(dropped_mis_barcode_table,
                                           as.numeric))
  rownames(mb_table_numeric) <- dropped_row_names
  colnames(mb_table_numeric) <- dropped_col_names
  #taxonomy now
  message("Starting taxonomic table building (long...)")
  mb.taxonIDneeded <- as.data.frame(row.names(mb_table_numeric))
  setnames(mb.taxonIDneeded, "row.names(mb_table_numeric)", "taxID")
  mb.taxonIDneeded <- as.numeric(mb.taxonIDneeded$taxID)
  message("Now downloading and putting together the NCBI
          database. This might take a while...")
  prepareDatabase(getAccessions=FALSE, indexTaxa=TRUE)
  message("Assigning all taxID in count matrix to fill taxonomy,
          you might want to take a break.")
  full.taxon.wimp <- getTaxonomy(mb.taxonIDneeded,'nameNode.sqlite')
  full.taxon.wimp.dt <- as.data.table(full.taxon.wimp, keep.rownames = "taxID")
  full.taxon.wimp.dt$taxID <- as.numeric(full.taxon.wimp.dt$taxID)
  full.taxon.wimp.dt$taxID <- as.numeric(full.taxon.wimp.dt$taxID)
  mb.dt <- as.data.table(mb_table_numeric, keep.rownames = "taxID")
  mb.dt$taxID <- as.numeric(mb.dt$taxID)
  merged.wimp.data <- merge(x = mb.dt, y = full.taxon.wimp.dt,
                            by = "taxID", all.x = TRUE)
  message("Done with taxa assignments. It should not be too much longer...")

  #get rid of unclassified
  if (keep.unclassified) {
    wimp.data.unclass.flag <- merged.wimp.data
  } else {
    wimp.data.unclass.flag <- merged.wimp.data[phylum != "NA"]
  }

  #get rid of human

  if (keep.human) {
    wimp.data.filtered <- wimp.data.unclass.flag
  } else {
    wimp.data.filtered <- wimp.data.unclass.flag[taxID != "9606"]
  }

  #prepare for phyloseq
  taxa_long <- wimp.data.filtered[, c("taxID", "superkingdom", "phylum", "class",
                                      "order", "family", "genus", "species")]
  count_table_mb <- wimp.data.filtered[, !c("superkingdom", "phylum", "class",
                                            "order", "family", "genus",
                                            "species")]
  mb.metadata <- as.data.frame(mb.metadata)
  rownames(mb.metadata) <- mb.metadata$sampleID

  OTU = otu_table(count_table_mb, taxa_are_rows = TRUE)
  TAX = tax_table(as.matrix(taxa_long))
  META = sample_data(mb.metadata)
  phyloseq(OTU, TAX, META)
}
