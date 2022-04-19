#' epi2me2r.
#'
#' @name epi2me2r
#' @description data for CARD db
#' @docType package
NULL

#' CARD 3.1.4 and 1.1.3 aro_index hand curated and combine
#'
#' A dataset the information on CARD genes in the 3.1.4 and 1.13 data releases
#' The variables are as follows:
#'
#' \itemize{
#'   \item ARO Accession. ARO Accession (3000005--3005385)
#'   \item CARDversion. If the data came from 3.1.4 or 1.1.3
#'   (3.1.3, hand-added1.1.3)
#'   \item CVTERMID CV TERM ID (36014--43745)
#'   \item Model Sequence ID. Model Sequence ID
#'   \item Model ID.
#'   \item Model Name.
#'   \item ARO Name.
#'   \item Protein Accession.
#'   \item DNA Accession.
#'   \item AMR Gene Family.
#'   \item Drug Class.
#'   \item Resistance Mechanism.
#'   \item mutation-associated. if the gene is mutation associated
#' }
#'
#' @docType data
#' @keywords datasets
#' @name CARD_taxonomy
#' @usage data(CARD_taxonomy)
#' @format A data frame with 3385 rows and 13 variables
#' @source \url{http://card.mcmaster.ca/download}
NULL
