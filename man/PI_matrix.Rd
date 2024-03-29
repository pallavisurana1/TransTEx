% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PI_calculation.R
\name{PI_matrix}
\alias{PI_matrix}
\title{Calculate Proportion Index (PI) Matrix}
\usage{
PI_matrix(data, iter = 0.5)
}
\arguments{
\item{data}{A data frame containing sample and ENST columns with values to be tested against iter.}

\item{iter}{A numeric threshold for determining whether the value in ENST columns is counted (default is 0.5).}
}
\value{
A data frame with proportions for each tissue and transcript (ENST columns).
}
\description{
This function calculates the Proportion Index matrix by first determining the count of values greater
than or equal to a threshold (iter) for each transcript (ENST columns) grouped by sample class,
then merging with tissue data and converting counts to proportions.
}
