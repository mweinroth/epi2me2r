#' Generates a taxonomy file for WIMP court table using taxonomizr
#' @name generate_wimp_taxonomy
#' @param wimp.count.table count of WIMP genes with taxID in the
#' first row and samples on columns
#' @return data table of taxonomy of WIMP genes in count table
#' @examples
#' \dontrun{
#' generate_wimp_taxonomy(wimp.count.table = wimp.count.table)
#' }
#' @import taxonomizr
#' @import data.table
#' @export
#'

generate_wimp_taxonomy <- function(wimp.count.table){
  # change to lower and always use lowercase to get columns
  mb.taxonIDneeded <- as.numeric(wimp.count.table$taxid)
  message("Now downloading and putting together the NCBI database. This might take a while...")
  prepareDatabase(getAccessions=FALSE, indexTaxa=TRUE)
  message("Assigning all taxID in count matrix to fill taxonomy. You might want to take a break...")
  full.taxon.wimp <- getTaxonomy(mb.taxonIDneeded,'nameNode.sqlite')
  full.taxon.wimp.dt <- as.data.table(full.taxon.wimp, keep.rownames = "taxid")
  full.taxon.wimp.dt$taxid <- as.numeric(full.taxon.wimp.dt$taxid)
  wimp.count.table$taxid <- as.numeric(wimp.count.table$taxid)
  merged.wimp.data <- merge(x = wimp.count.table, y = full.taxon.wimp.dt,
                            by = "taxid", all.x = TRUE)
  taxa_long <- merged.wimp.data[, c("taxid", "superkingdom", "phylum",
                                    "class", "order", "family", "genus",
                                    "species")]
  taxa_long
}
