\name{CLX.test.unequalcov}
\alias{CLX.test.unequalcov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cai, Liu, and Xia equal means test under unequal covariances
%%  ~~function to do ... ~~
}
\description{
	
	Performs the test in Cai, Liu, and Xia (2014) for the equality of two \code{p} by \code{1} population mean vectors given samples of sizes \code{n} and \code{m} when the popoulation covariance matrices are assumed to be unequal.	
	
	}
	
\usage{
CLX.test.unequalcov(X, Y)
}

\arguments{
  \item{X}{ the \code{n} by \code{p} data matrix for sample one.}
  \item{Y}{ the \code{m} by \code{p} data matrix for sample two.}
}
\value{
\item{TSvalue}{the value of the test statistic.}
\item{pvalue}{the two-sided p-value for the test statistic.}
}
\references{ 

	Cai, T. T., Liu, W. & Xia, Y. (2014). Two-sample test of high dimensional means under dependence.
		\emph{J. R. Statist. Soc. B.}
}
\author{
Karl Gregory   \email{kgregory@mail.uni-mannheim.de}, \url{http://www.stat.tamu.edu/~kbgregory}.
}

\seealso{
\code{\link{CLX.sim.equalcov}}
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
	
	CLX.test.unequalcov(X,Y)
	
}
}
\keyword{htest}
\keyword{multivariate} 
