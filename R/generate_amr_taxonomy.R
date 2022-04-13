#' Generates a taxonomy file for AMR court table
#' @name generate_amr_taxonomy
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

generate_amr_taxonomy <- function(amr.count.table, verbose=FALSE){
  amr.CVTERMID.list <- amr.count.table[,1]
  amr.CVTERMID.list$CVTERMID <- as.numeric(amr.CVTERMID.list$CVTERMID)
  CARD_taxonomy$CVTERMID <- as.numeric(CARD_taxonomy$CVTERMID)
  merged.data <- merge(x = amr.CVTERMID.list, y = CARD_taxonomy, by = "CVTERMID", all.x = TRUE)
  if (verbose == TRUE) {
    taxa_long <- merged.data[, c("CVTERMID", "ARO Accession", "CARDversion", "Model Sequence ID", "Model ID",
                                 "Model Name", "ARO Name", "Protein Accession", "DNA Accession",
                                 "AMR Gene Family", "Drug Class", "Resistance Mechanism", "mutationassociated")]
    taxa_long}
  else if(verbose == FALSE) {
    taxa_short <- merged.data[, c("CVTERMID","Drug Class", "AMR Gene Family", "Resistance Mechanism", "ARO Name")]
    taxa_short}
}
