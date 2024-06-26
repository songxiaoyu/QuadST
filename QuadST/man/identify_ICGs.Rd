% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QuadST-main-functions.R
\name{identify_ICGs}
\alias{identify_ICGs}
\title{Identify cell-cell interaction changes genes (IGGs)}
\usage{
identify_ICGs(
  x,
  y,
  dist,
  expr,
  cov = NULL,
  tau,
  p_thres = 0.05,
  fdr = 0.1,
  ABconst = 0.1
)
}
\arguments{
\item{x}{A \code{SingleCellExperiment} class.}

\item{y}{A matrix of quntile regression p-values: genes (rows) by quantiles (columns).}

\item{dist}{A column name of \code{colData(object1)} that stores anchor-neighbor cell pair distances.}

\item{expr}{A column name of \code{assays(object1)} that stores anchor cells' gene expression levels.}

\item{cov}{Column names of \code{colData(object1)} that need to be adjusted as covariates.}

\item{tau}{A set of quantile levels at which test statistics are calculated.}

\item{p_thres}{An initial p-value threshold value. Use 0.05 by default. This serves as a starting point to find a significant p-value with the empirical FDR procedure.}

\item{fdr}{A norminal false discovery rate to control. Use 0.1 by default.}

\item{ABconst}{A constant used in the empirical FDR procedure. Use 0.1 by default.}
}
\value{
a list of inference results
\item{ACAT}{Combined gene-specific p-values at each lowest and highest quantile symmetrically around the median quantile using Cauchy Combination test.}
\item{p_sig}{Significant p value threshold across quantile levels at a prescribed FDR threshold, e.g., 0.1.}
\item{q_sig}{Logicals that indicate if there exist significant p value thresholds across quantile levels.}
\item{ICGs}{Identified cell-cell interaction changed gene names.}
\item{q_int}{An interaction quantile.}
\item{dist_int}{An interaction distance.}
\item{DA_score}{Directional association scores calculated at an interaction quantile.}
}
\description{
Identify cell-cell interaction changes genes (IGGs)
}
