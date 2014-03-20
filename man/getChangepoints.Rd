\name{getChangepoints}
\alias{getChangepoints}
\title{Getting Change-points}
\description{This function estimates the change-points where one or both parent-specific copy numbers change.  It uses a circular binary segmentation approach to find change-points in a binomial mixture process.  The output of the function is the set of locations of the break points.  If the whole genome is analyzed, it is recommended to run this function chromosome by chromosome, and the runs on different chromosomes can be done in parallel to shorten the running time.
}
\usage{
getChangepoints(AT, BT, AN, BN, verbose=TRUE, pOri=c(0.49,0.51), error=1e-5, maxIter=1000)
}
\arguments{
  \item{AT}{A vector of A-allele coverage at the variant loci in the case sample.}
  \item{BT}{A vector of B-allele coverage at the variant loci in the case sample.}
  \item{AN}{A vector of A-allele coverage at the variant loci in the control sample.}
  \item{BN}{A vector of B-allele coverage at the variant loci in the control sample.}
  \item{verbose}{Provide progress messages if it is TRUE.  This argument is TRUE by default.  Set it to be FALSE if you want to turn off the progress messages.}
  \item{pOri, error, maxIter}{Parameters used in estimating the success probabilities of the mixed binomial distribution.  See the manuscript by Chen and Zhang for more details.  "pOri" provides the initial success probabilities.  The two values in pOri needs to be different.  "error" provides the stopping criterion.  "maxIter" is the maximum iterating steps if the stopping criterion is not achieved.}
}
\seealso{
  \code{\link{getASCN}}
}
\references{
Chen, H. and Zhang, N.R. (2014) Allele-specific copy number profiling by next-generation DNA sequencing.
}
\examples{
data(Example) 
attach(ReadCounts)
# tauhat = getChangepoints(AT, BT, AN, BN)  # uncomment this to run the function. 
 # This example has 6309 variant loci and it took 94 seconds to run on a laptop with 
 # Intel Core i5-2410M processor. 
}

