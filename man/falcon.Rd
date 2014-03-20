\name{falcon}
\alias{falcon}
\title{Finding Allele-specific Copy Number in Next-Generation Sequencing Data}
\description{This library contains a set of tools for Allele-specific DNA Copy Number Profiling using Next-Generation Sequencing.  Given the allele-specific coverage at a set of variant loci, get getChangepoints program segments the genome into regions of homogeneous allele-specific copy number.  It requires, as input, the read counts supporting each variant allele in each of a pair of case and control samples.  For example, for detection of somatic mutations in a tumors, the case and control samples can be, respectively, the tumor and normal sample from the same individual.
}
\author{
Hao Chen and Nancy R. Zhang

Maintainer: Hao Chen (hxchen@ucdavis.edu)
}
\references{
Chen, H. and Zhang, N.R. (2014) Allele-specific copy number profiling by next-generation DNA sequencing.
}

\seealso{
  \code{\link{getChangepoints}}, \code{\link{getASCN}}, \code{\link{view}}
}
\examples{
data(Example) 
attach(ReadCounts)
# tauhat = getChangepoints(AT, BT, AN, BN)  # uncomment this to run the function. 
 # This example has 6309 variant loci and it took 94 seconds to run on a laptop with 
 # Intel Core i5-2410M processor.
cn = getASCN(AT, BT, AN, BN, tauhat=tauhat)
 # cn$tauhat would give the indices of change-points.  
 # cn$ascn would give the estimated allele-specific copy numbers for each segment.
 # cn$Haplotype[[i]] would give the estimated haplotype for the major chromosome in segment i 
 # if this segment has different copy numbers on the two homologous chromosomes.
view(AT,BT,AN,BN,cn$tauhat,cn$ascn)
}

