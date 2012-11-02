\name{SK.test}
\alias{SK.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Srivastava and Kubokawa test for equal means
%%  ~~function to do ... ~~
}
\description{
	
	Performs the test from Srivastava and Kubokawa (2013) for the equality of two \code{p} by \code{1} population mean vectors given samples of sizes \code{n} and \code{m}.	
	
	}
	
\usage{
SK.test(X, Y)
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

	Srivastava, M. S. and Kubokawa, T. (2013). Tests for multivariate analysis of variance in high dimension under non-normality. 
		\emph{Journal of Multivariate Analysis} \bold{115}, 204216.	
}
\author{
Karl Gregory   \email{kgregory@mail.uni-mannheim.de}, \url{http://www.stat.tamu.edu/~kbgregory}.
}

\seealso{
\code{\link{SK.sim}}
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
	
	SK.test(X,Y)
	
	}
}
\keyword{htest}
\keyword{multivariate} 
