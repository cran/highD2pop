\name{chr1qseg}
\alias{chr1qseg}
\title{Copy number data along a segment of the q arm of chromosome 1}
\description{.Rdata object containing copy number probes from TCGA at 400 locations (given in 'loc', which is number of base pairs from end of chromosome arm) on the q arm of chromosome 1 for 92 long term survivors 'X' and 138 short term survivors 'Y' of glioblastoma multiforme.  Around 3\% of the values are missing.}
\usage{
chr1qseg
}

\format{ A list containing \code{X}, the \code{n=92} by \code{p=400} matrix of copy number measurements for the \code{n=92} long-term survivors, \code{Y}, the \code{m=138} by \code{p=400} matrix of copy number measurements for the \code{m=138} short-term survivors, and \code{loc}, which gives the locations of the \code{400} measurements in number of base pairs from the end of the q arm of chromosome 1.
	}	


\keyword{datasets}
