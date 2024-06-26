% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QuadST-main-functions.R
\name{test_QuadST_model}
\alias{test_QuadST_model}
\title{Test anchor-neighbor distance-expression association
 at a set of highest and lowest quantiles symmetric around median}
\usage{
test_QuadST_model(x, dist, expr, cov = NULL, tau)
}
\arguments{
\item{x}{A \code{SingleCellExperiment} class.}

\item{dist}{A column name of \code{colData(object)} that stores anchor-neighbor cell pair distances.}

\item{expr}{A column name of \code{assays(object)} that stores anchor cells' gene expression levels.}

\item{cov}{Column names of \code{colData(object)} that needs to be adjusted as covariates.}

\item{tau}{A set of highest and lowest quantiles symmetric around median.}
}
\value{
A matrix of quntile regression p-values: genes (rows) by quantiles (columns).
}
\description{
Test anchor-neighbor distance-expression association
 at a set of highest and lowest quantiles symmetric around median
}
