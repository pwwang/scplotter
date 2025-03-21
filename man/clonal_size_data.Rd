% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clonalutils.R
\name{clonal_size_data}
\alias{clonal_size_data}
\title{clonal_size_data}
\usage{
clonal_size_data(data, clone_call, chain, groupings)
}
\arguments{
\item{data}{The product of \link[scRepertoire:combineTCR]{scRepertoire::combineTCR}, \link[scRepertoire:combineTCR]{scRepertoire::combineTCR}, or
\link[scRepertoire:combineExpression]{scRepertoire::combineExpression}.}

\item{clone_call}{How to call the clone - VDJC gene (gene), CDR3 nucleotide (nt),
CDR3 amino acid (aa), VDJC gene + CDR3 nucleotide (strict) or a custom variable
in the data}

\item{chain}{indicate if both or a specific chain should be used - e.g. "both",
"TRA", "TRG", "IGH", "IGL"}

\item{groupings}{The column names in the meta data to group the cells.}
}
\value{
A data frame with the clonal size data.
}
\description{
Function to get the clonal size data for all group_by values.
}
\keyword{internal}
