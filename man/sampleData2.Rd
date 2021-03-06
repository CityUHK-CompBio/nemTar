% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampleData2.R
\name{sampleData2}
\alias{sampleData2}
\title{Generation of the data for sampling in the simulation of the cancer multi-omics samples}
\usage{
sampleData2(
  Phi,
  S_obs,
  p,
  m,
  type = "binary",
  typeI.err = typeI.err,
  typeII.err = typeII.err
)
}
\arguments{
\item{Phi}{adjacency matrix}

\item{S_obs}{binary matrix inddicating the observation of S-genes' states}

\item{p}{number of samples(patients) to sample}

\item{m}{number of E-genes to sample}

\item{type}{an integer,represents the number of row within the raw image}

\item{typeI.err}{simulated type I error for binary data}

\item{typeII.err}{simulated type II error for binary data}
}
\description{
Generation of the data for sampling in the simulation of the cancer multi-omics samples
}
\seealso{
\code{\link[nem]{sampleData}}
}
