##' Raw WIMP files plus metadata to metagenomeseq object
##'@name wimp_raw_to_metagenomeseq
##' @description Given WIMP directory and metadata, make phyloseq
##' object. \pkg{\link{metagenomeSeq}} package required.
##' @param path.to.wimp.files path to data of raw CSV files from
##' WIMP analysis
##' @param metadata dataframe of metadata with "filename" and "barcode"
##' columns required
#' @param keep.unclassified TRUE or FALSE: whether to keep reads that do not
#' classify below phylum, default = FALSE
#' @param keep.human TRUE or FALSE: whether to keep reads that are classified as
#' human, default = FALSE
##' @seealso \pkg{\link{metagenomeSeq}}
#' @return metagenomeseq object for downstream analysis with WIMP data
#' @examples
#' \dontrun{
#' wimp_raw_to_metagenomeseq(path.to.wimp.files = path/to/wimpfiles,
#' metadata = metadata, keep.unclassified = FALSE, keep.human = FALSE)
#' }
#' @import data.table
#' @import taxonomizr
#' @importFrom  phyloseq phyloseq
#' @importFrom  phyloseq otu_table
#' @importFrom phyloseq tax_table
#' @importFrom  phyloseq sample_data
#' @importFrom Biobase AnnotatedDataFrame
#' @import metagenomeSeq
#' @export


wimp_raw_to_metagenomeseq <- function(path.to.wimp.files, metadata,
                                      keep.unclassified = FALSE,
                                      keep.human = FALSE){

  # Checks for valid input. Fails with error if any are not met.
  stopifnot(dir.exists(path.to.amr.files))
  stopifnot(is.data.frame(metadata))
  stopifnot(is.logical(keep.unclassified))
  stopifnot(is.logical(keep.human))

  if (any(!c('filename', 'barcode') %in% names(metadata))) {
    stop('Error: metadata does not have columns named "filename" and "barcode".')
  }

  #read in raw files
  message(paste("Reading in raw files from", path.to.wimp.files))
  parsed_files <- list.files(path = path.to.wimp.files)
  Sample_IDs <- sub(".csv", "", parsed_files)

  mb.rawdata <- lapply(1:length(parsed_files), function(i) {
    file_name <- paste0(parsed_files[i])
    mb.dataframe <- fread(paste0(path.to.wimp.files,file_name))
    cbind(mb.dataframe, csvname = Sample_IDs[i])
  })

  mb.rawdata <- do.call(rbind, mb.rawdata)

  total.reads <- nrow(mb.rawdata)
  mb.classified <- mb.rawdata[mb.rawdata$exit_status == "Classified",]
  classified.reads <- nrow(mb.classified)
  percentage.classified <- round((classified.reads/total.reads*100),
                                digits = 2)
  message(paste("The percentage of classified reads was",
                percentage.classified, "%"))
  mb.rawdata.reduced <- mb.rawdata[, list(csvname, barcode, taxID)]
  sampleidinfo <- mb.rawdata.reduced[, .(sampleID=
                                           mb.rawdata.reduced[["sampleID"]],
                                         sampleID=
                                           do.call(paste, c(.SD, sep="_"))),
                                     .SDcols= csvname:barcode]
  mb.rawdata.reduced <- cbind(sampleidinfo, mb.rawdata.reduced)
  mb.rawdata.reduced.subset <- mb.rawdata.reduced[, list(sampleID, taxID)]
  wimp_count_table <- suppressMessages(dcast(mb.rawdata.reduced.subset,
                                             taxID ~ sampleID))
  #remove barcodes not in metadata
  #metadata[] <- lapply(metadata, function(y) gsub("BC", "barcode", y))
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
  #entire section because I could not get the df to be numeric in the count data
  dropped_mis_barcode <- as.data.frame(t(dropped_mis_barcode.t))
  dropped_col_names <- as.character(dropped_mis_barcode[1,])
  dropped_mis_barcode_table <- dropped_mis_barcode[-1, ]
  dropped_row_names <- as.character(row.names(dropped_mis_barcode_table))
  mb_table_numeric <- as.data.frame(sapply(dropped_mis_barcode_table,
                                           as.numeric))
  rownames(mb_table_numeric) <- dropped_row_names
  colnames(mb_table_numeric) <- dropped_col_names
  #taxonomy now
  message("starting taxonomic table building (long...)")
  mb.taxonIDneeded <- as.data.frame(row.names(mb_table_numeric))
  setnames(mb.taxonIDneeded, "row.names(mb_table_numeric)", "taxID")
  mb.taxonIDneeded <- as.numeric(mb.taxonIDneeded$taxID)
  message("Now downloading and putting together the NCBI database.
          This might take a while...")
  prepareDatabase(getAccessions=FALSE, indexTaxa=TRUE)
  message("Assigning all taxID in count matrix to fill taxonomy.
          You might want to take a break...")
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
  if (keep.unclassified == TRUE) {
    wimp.data.unclass.flag <- merged.wimp.data
  } else if (keep.unclassified == FALSE) {
    wimp.data.unclass.flag <- merged.wimp.data[phylum != "NA"]
  }

  #get rid of human
  if (keep.human == TRUE) {
    wimp.data.filtered <- wimp.data.unclass.flag
  } else if (keep.human == FALSE) {
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
  ps <- (phyloseq(OTU, TAX, META))
  phyloseq_to_metagenomeSeq(ps)
}
