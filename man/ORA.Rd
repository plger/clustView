% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plg.GO.R
\name{ORA}
\alias{ORA}
\title{ORA}
\usage{
ORA(set1, sets, universe, topN = 10)
}
\arguments{
\item{set1}{Character vector; the set to test}

\item{sets}{Named list containing the sets for which to test enrichment}

\item{universe}{A character vector containing the universe/background, 
or a numeric or length 1 containing the size of the universe.}

\item{topN}{The number of top enrichments to return (default 10)}
}
\value{
A data.frame
}
\description{
A simple hypergeometric over-enrichment analysis
}
