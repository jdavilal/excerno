% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classify_simluated_samples.R
\name{classify_simulated_samples}
\alias{classify_simulated_samples}
\title{Get classification of each mutation}
\usage{
classify_simulated_samples(samples, signatures, signature.names = c())
}
\arguments{
\item{samples}{List of mutation samples, each coming from a different mutational signature}

\item{signatures}{List of the mutational signature that each sample came from. Each element should correspond an element in argument samples at the same index.}

\item{signature.names}{Vector of the signature names. Each element should correspond an element in argument signatures at the same index.}
}
\value{
A data frame with the classification of each mutation in the samples
}
\description{
Get classification of each mutation using Bayes' Theorem.
}
\examples{

library(MutationalPatterns)
library(tidyverse)

# Load in signatures
cosmic.sigs <- get_known_signatures()
cosmic.sig4 <- as.matrix(cosmic.sigs[,4])
ffpe.sig <- get_ffpe_signature()

# Create samples
sample.sig4 <- create_signature_sample_vector(cosmic.sig4, 100)
sample.ffpe <- create_signature_sample_vector(ffpe.sig, 100)

samples <- list(sample.sig4, sample.ffpe)
signatures <- list(cosmic.sig4, ffpe.sig)

classify_simulated_samples(samples, signatures)
classify_simulated_samples(samples, signatures, c("SBS4", "FFPE"))
}
