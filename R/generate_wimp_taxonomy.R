#' Generates a taxonomy file for WIMP court table using taxonomizr
#' @name generate_wimp_taxonomy
#' @param wimp.count.table count of WIMP genes with taxID in the first row adn samples on columns
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
  mb.taxonIDneeded <- as.numeric(wimp.count.table$taxID)
  message("Now downloading and putting together the NCBI databasem this might take a while...")
  prepareDatabase(getAccessions=FALSE, indexTaxa=TRUE)
  message("Assigning all taxID in count matix to fill taxonomy, you might want to take a break.")
  full.taxon.wimp <- getTaxonomy(mb.taxonIDneeded,'nameNode.sqlite')
  full.taxon.wimp.dt <- as.data.table(full.taxon.wimp, keep.rownames = "taxID")
  full.taxon.wimp.dt$taxID <- as.numeric(full.taxon.wimp.dt$taxID)
  wimp.count.table$taxID <- as.numeric(wimp.count.table$taxID)
  merged.wimp.data <- merge(x = wimp.count.table, y = full.taxon.wimp.dt, by = "taxID", all.x = TRUE)
  taxa_long <- merged.wimp.data[, c("taxID", "superkingdom", "phylum", "class", "order", "family", "genus", "species")]
  print(taxa_long)
}
