% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mpcross-class.R
\name{infiniteSelfing}
\alias{infiniteSelfing}
\title{Create allele encoding corresponding to infinite generations of selfing}
\usage{
infiniteSelfing(founders, finals, pedigree)
}
\arguments{
\item{founders}{The genetic data for the founding lines, which are assumed to be inbred}

\item{finals}{The genetic data for the lines genotyped at the end of the experiment.}

\item{pedigree}{The pedigree for the experiment}
}
\value{
An object of class \code{hetData}, which encodes only the marker homozygotes.
}
\description{
Create allele encoding corresponding to infinite generations of selfing
}
\details{
In many experiments (particularly those that are significantly inbred), only marker homozygotes are observed, which means that the relationship between marker genotypes and marker alleles is particularly simple. In such cases, generally a marker genotype of some value (say 0) indicates that the individual is homozygous for marker allele 0. 

This function takes in genetic data for the founding lines, genetic data for the final population, and the pedigree. It returns an encoding for marker genotypes where every genotype is homozygous for the marker allele with the same value.
}
\examples{
map <- qtl::sim.map()
pedigree <- f2Pedigree(1000)
cross <- simulateMPCross(map = map, pedigree = pedigree, mapFunction = haldane, seed = 1)
#Initially the object contains markers that are fully informative.
#The final genetic data contains values 1, 2 and 3, while the genetic data for the founding 
#    lines contains only values 1 and 2. 
#A value of 1 or 2 in the final genetic data indicates a homozygote for the 
#	corresponding marker allele. 
#A value of 3 in the final genetic data indicates a heterozygote for the marker allele.
#Information about this encoding is stored in the hetData slot.
hetData(cross, "D1M1")
cross <- cross + biparentalDominant()
#Now we have converted all markers to dominant.
#The final genetic data contains values 1 and 2, and the genetic data for the founding 
#    lines contains only values 1 and 2. 
#A value of 2 indicates a homozygote for the corresponding marker allele, OR a 
#	marker heterozygote.
hetData(cross, "D1M1")
#But under infinite generations of selfing, the encoding is simpler. 
simpleEncoding <- infiniteSelfing(founders = founders(cross), finals = finals(cross), 
	pedigree = pedigree)
simpleEncoding[["D1M1"]]
}
