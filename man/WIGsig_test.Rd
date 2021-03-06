% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WIGsig_test.R
\name{WIGsig_test}
\alias{WIGsig_test}
\title{The statistical significance test based on the bootstrap of WIGs}
\usage{
WIGsig_test(sample_WIG0, WIG, sampling_times)
}
\arguments{
\item{sample_WIG0}{the background distribution of the WIG from the sampling of WIG}

\item{WIG}{the pathway-specific WIGcalculated using \code{compute_WIG}}

\item{sampling_times}{a numeric number indicating the times of sampling}
}
\value{
the statistical significance(BH-adjusted p-value)
}
\description{
The statistical significance test based on the bootstrap of WIGs
}
