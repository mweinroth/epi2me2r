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

  ps <- amr_raw_to_phyloseq(path.to.amr.files, metadata, coveragenumber, keepSNP)
  phyloseq_to_metagenomeSeq(ps)
}
