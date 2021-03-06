% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getIntercrossingAndSelfingGenerations.R
\name{getIntercrossingAndSelfingGenerations}
\alias{getIntercrossingAndSelfingGenerations}
\title{Identify number of generations of intercrossing and selfing, per genetic line}
\usage{
getIntercrossingAndSelfingGenerations(cross)
}
\arguments{
\item{cross}{The \code{mpcross} object containing the pedigree to be analysed.}
}
\value{
An integer matrix with two columns, giving the number of generations of selfing and intercrossing, for each genetic line. Or in the case of multiple experiments contained within a single object, a list of such matrices.
}
\description{
Identify number of generations of intercrossing and selfing, per genetic line
}
\details{
Many structured populations consist of a number of generations of mixing, followed by a number of generations of intercrossing, followed by inbreeding. This function identifies the number of generations of selfing and intercrossing, for each genetic line, in the case of 4-way, 8-way or 16-way multi-parent design.
}
