% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p_screen.R
\name{p_screen}
\alias{p_screen}
\title{Compute screening P-value for one gene}
\usage{
p_screen(Z, Sigma)
}
\arguments{
\item{Z}{numeric, vector of TWAS Z-scores}

\item{Sigma}{numeric, matrix of TWAS LD within a gene}
}
\value{
screen P-value for one gene
}
\description{
The function runs a modified version of P_ACT to test
hypothesis of at least one gene-transcript association
}
