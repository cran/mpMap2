% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/impute.R
\name{impute}
\alias{impute}
\title{Impute missing recombination fraction estimates}
\usage{
impute(
  mpcrossLG,
  verbose = FALSE,
  allErrors = FALSE,
  extractErrorsFunction = function(e) e
)
}
\arguments{
\item{mpcrossLG}{An object of class \code{mpcrossLG}, which contains estimated pairwise recombination fractions}

\item{verbose}{Should more verbose output be generated?}

\item{allErrors}{If there is an error, should we immediately return, or should we continue, and report all errors?}

\item{extractErrorsFunction}{Error handling function. If there are errors and allErrors is \code{TRUE}, this function will be called with a matrix indicating which estimates could not be imputed.}
}
\value{
An object of class \code{mpcrossLG}, containing all the information in the input object, but also an imputed copy of the estimated recombination fraction data.
}
\description{
Impute missing recombination fraction estimates
}
\details{
Recombination fractions between every pair of markers are estimated using numerical maximum likelihood. Unfortunately the likelihood is flat in some cases, so an estimate cannot be made. This later causes problems when trying to use estimated recombination fractions to order the markers, because a complete matrix of estimates is required. The solution is to impute the missing estimates using related estimates. For example, the recombination fraction between markers A and C may not be directly estimatable. However, there may be a marker B known to be tightly linked to A, which has a known recombination fraction with C. The estimated recombination fraction between B and C can be taken to be an estimate of the recombination fraction between A and C.

This function imputes values in the estimated recombination fraction matrix, to return a complete matrix. If there is a value that cannot be imputed, an error is triggered. Input \code{allErrors} controls whether the function will stop after encountering a single error, or continue and report all errors. If all errors are being reported, the optional function \code{extractErrorsFunction} is called with information about which missing estimates could not be imputed.
}
