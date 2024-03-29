% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wimp_raw_to_phyloseq.R
\name{wimp_raw_to_phyloseq}
\alias{wimp_raw_to_phyloseq}
\title{Raw WIMP files plus metadata to phyloseq object}
\usage{
wimp_raw_to_phyloseq(
  path.to.wimp.files,
  metadata,
  keep.unclassified = FALSE,
  keep.human = FALSE
)
}
\arguments{
\item{path.to.wimp.files}{path to data of raw csv files from WIMP
analysis}

\item{metadata}{dataframe of metadata with "filename" and "barcode"
columns required}

\item{keep.unclassified}{TRUE or FALSE: whether to keep reads that do not
classify below phylum, default = FALSE}

\item{keep.human}{TRUE or FALSE: whether to keep reads that are classified as
human, default = FALSE}
}
\value{
phyloseq object for downstream analysis with WIMP data
}
\description{
Given wimp directory and metadata, make phyloseq object.
\pkg{\link{phyloseq}} package required.
}
\examples{
\dontrun{
wimp_raw_to_phyloseq(path.to.wimp.files = path/to/wimpfiles,
metadata = metadata, keep.unclassified = FALSE, keep.human = FALSE)
}
}
\seealso{
\pkg{\link{phyloseq}}
}
