\name{CLX.Covtest}
\alias{CLX.Covtest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cai, Liu, and Xia test for equality of two mean vectors
%%  ~~function to do ... ~~
}
\description{
	
	Performs the test from Cai, Liu, and Xia (2013) for the equality of two population covariance matrices.
	}
	
\usage{
CLX.Covtest(X, Y)
}

\arguments{
  \item{X}{ the \code{n} by \code{p} data matrix for sample one.}
  \item{Y}{ the \code{m} by \code{p} data matrix for sample two.}
}
\value{
\item{TSvalue}{the value of the test statistic.}
\item{pvalue}{the p-value for the test statistic.}
}
\references{ 
	T. Tony Cai, Weidong Liu, and Yin Xia. 2013 
	``Two-sample covariance matrix testing and support recovery in high-dimensional and sparse settings''  		
	\emph{Journal of The American Statistical Association}. \bold{108}, 265--277
}
\author{
Karl Gregory   \email{kgregory@mail.uni-mannheim.de}, \url{http://www.stat.tamu.edu/~kbgregory}.
}

\seealso{
\code{\link{CLX.sim.Covtest}}
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
	
	CLX.test.equalcov(X,Y)
	
	}
}
\keyword{htest}
\keyword{multivariate} 
