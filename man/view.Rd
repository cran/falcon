\name{view} \alias{view} \title{Viewing Data with Allele-specific Copy Number} 

\description{This function generates three plots:  The first plots the A-allele frequencies of the case (black) sample overlayed onto those of the control (gray) sample; the second plots the relative depth of the case over control adjusted by the ratio of total mapped 
reads, i.e. P*(read count in tumor)/(read count in normal), where P=(total reads mapped in normal)/(total reads mapped in tumor); the third plots the estimated 
parent-specific DNA copy numbers. } 
\usage{ view(AT, BT, AN, BN, tauhat, ascn, pos=NULL, rdep=NULL, plot="all", ...) } 
\arguments{ 
  \item{AT}{A vector of A-allele reads in the case sample.}
  \item{BT}{A vector of B-allele reads in the case sample.}
  \item{AN}{A vector of A-allele reads in the control sample.}
  \item{BN}{A vector of B-allele reads in the control sample.}
  \item{tauhat}{A vector containing the estimated break points.  If it is not specified (NULL), then this function will first estimate the break points by calling the function "gettau", and then estimate the parent-specific DNA copy number for each segment.}
  \item{ascn}{The estimated parent-specific DNA copy numbers.}
  \item{pos}{A vector of the base positions for the SNPs.  If this information is not provided, the x-axis of the plots will simply be the SNP ordering.  If this information is provided, the x-axis of the plots will be the position information.} 
  \item{rdep}{The relative depth of the case sample over the control sample.  If it is not specified (NULL), then the value median(AT+BT)/median(AN+BN) will be used.}
  \item{plot}{This argument determines what to plot.  By default, this function gives all three plots described above ("all").  You can also plot each one individually if you set this argument to either of "Afreq", "RelativeCoverage" or "ASCN".}
  \item{...}{Arguments from plot can be passed along.}
}
\seealso{
  \code{\link{getASCN}}
}
\examples{
data(Example) 
attach(ReadCounts)
cn = getASCN(AT, BT, AN, BN, tauhat=tauhat)
view(AT,BT,AN,BN,cn$tauhat,cn$ascn)

# to view the plot for only showing A-allele frequency of the case (black) sample overlayed 
# onto those of the control (gray) sample
par(mfrow=c(1,1))
view(AT,BT,AN,BN,cn$tauhat,cn$ascn, plot="Afreq")

# to view the relative depth of the case over control adjusted by the ratio of total mapped 
# reads in fixed size bins
par(mfrow=c(1,1))
view(AT,BT,AN,BN,cn$tauhat,cn$ascn, plot="RelativeCoverage")

# to view the estimated allele-specific DNA copy numbers
par(mfrow=c(1,1))
view(AT,BT,AN,BN,cn$tauhat,cn$ascn, plot="ASCN")
}

