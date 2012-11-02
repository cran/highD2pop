\name{ChenQin.test}
\alias{ChenQin.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Chen's and Qin's test for equality of two mean vectors
%%  ~~function to do ... ~~
}
\description{
	
	Performs the test from Chen and Qin (2010) for the equality of two \code{p} by \code{1} population mean vectors given samples of sizes \code{n} and \code{m}.	
	
	}
	
\usage{
ChenQin.test(X, Y)
}

\arguments{
  \item{X}{ the \code{n} by \code{p} data matrix for sample one.}
  \item{Y}{ the \code{m} by \code{p} data matrix for sample two.}
}
\value{
\item{ChQ}{the value of the test statistic.}
\item{pvalue}{the two-sided p-value for the test statistic.}
}
\references{ 
	Chen, X. S. and Qin, Y.L. (2010). A two sample test for high dimensional data with applications to gene-set testing.  		
		\emph{The Annals of Statistics}. \bold{38}(2):808--835
}
\author{
Karl Gregory   \email{kgregory@mail.uni-mannheim.de}, \url{http://www.stat.tamu.edu/~kbgregory}.
}

\examples{
\dontrun{
	
data(chr1qseg)
	
	impute <- function(x) 	
	{ 	
		x[which(is.na(x))] <- mean(x,na.rm=TRUE)
		return(x)
	}
	
	X <- apply(chr1qseg$X,2,impute)
	Y <- apply(chr1qseg$Y,2,impute)
	
	ChenQin.test(X,Y)
	
	}
}
\keyword{htest}
\keyword{multivariate} 
