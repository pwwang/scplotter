% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clonalutils.R
\name{merge_clonal_groupings}
\alias{merge_clonal_groupings}
\title{merge_clonal_groupings}
\usage{
merge_clonal_groupings(data, groupings, sep = " // ")
}
\arguments{
\item{data}{The product of \link[scRepertoire:combineTCR]{scRepertoire::combineTCR}, \link[scRepertoire:combineTCR]{scRepertoire::combineTCR}, or
\link[scRepertoire:combineExpression]{scRepertoire::combineExpression}.}

\item{groupings}{A list of the clonal groupings. Each element is a column in the data.}
}
\value{
The data with the combined groupings (\code{.group})
}
\description{
Merge the multiple clonal groupings into a single grouping.
}
\details{
Because \link[scRepertoire:clonalQuant]{scRepertoire::clonalQuant} and families don't support mutliple groupings,
this is trying to merge the multiple groupings into a single grouping. And then
later restore the original groupings.
}
\keyword{internal}
