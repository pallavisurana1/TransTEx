% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subsample.R
\name{calculate_empirical_p}
\alias{calculate_empirical_p}
\title{Calculate Empirical P-values}
\usage{
calculate_empirical_p(data, threshold = 0.05)
}
\arguments{
\item{data}{A data frame containing the dataset for analysis with 'tissue' and 'ENST' prefixed columns for transcripts.}

\item{threshold}{The p-value threshold used to determine significance (default at 0.05).}
}
\value{
A data frame with transcript IDs, empirical p-values, and other metadata if specified.
}
\description{
This function integrates the process of calculating sample sizes, performing stratified sampling and hypothesis testing,
and finally calculating empirical p-values for tissue-specific analysis.
}
