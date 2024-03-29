% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_wimp_taxonomy.R
\name{generate_wimp_taxonomy}
\alias{generate_wimp_taxonomy}
\title{Generates a taxonomy file for WIMP court table using taxonomizr}
\usage{
generate_wimp_taxonomy(wimp.count.table)
}
\arguments{
\item{wimp.count.table}{count of WIMP genes with taxID in the
first row and samples on columns}
}
\value{
data table of taxonomy of WIMP genes in count table
}
\description{
Generates a taxonomy file for WIMP court table using taxonomizr
}
\examples{
\dontrun{
generate_wimp_taxonomy(wimp.count.table = wimp.count.table)
}
}
