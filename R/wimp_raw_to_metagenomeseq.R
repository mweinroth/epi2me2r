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

  ps <- wimp_raw_to_phyloseq(path.to.wimp.files, metadata, keep.unclassified, keep.human)
  phyloseq_to_metagenomeSeq(ps)
}
