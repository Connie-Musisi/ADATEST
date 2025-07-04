% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Remove_NAs.R
\name{Remove_NAs}
\alias{Remove_NAs}
\title{Remove Taxa with Missing Values from a Phyloseq Object}
\usage{
Remove_NAs(physeq)
}
\arguments{
\item{physeq}{A phyloseq object containing OTU and taxonomy tables.}
}
\value{
A modified phyloseq object with taxa containing NAs removed.
}
\description{
This function removes taxa (OTUs) that contain missing values (NA) from a given phyloseq object.
It ensures that only taxa with complete data are retained in both the OTU and taxonomy tables.
}
