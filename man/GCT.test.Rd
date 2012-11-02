\name{GCT.test}
\alias{GCT.test}

\title{Generalized component test
}

\description{Performs the generalized component test from Gregory et al. (2014) for the equality of two \code{p} by \code{1} population mean vectors given samples of sizes \code{n} and \code{m}.
}

\usage{
GCT.test(X, Y, r, smoother = "parzen", ntoorderminus = 2)
}

\arguments{
  \item{X}{the \code{n} by \code{p} data matrix for sample one.
  	}
  \item{Y}{the \code{m} by \code{p} data matrix for sample two.
  	}
  \item{r}{the lag window size for variance estimation.
  	}
  \item{smoother}{the lag window used in the variance estimation.  Possible values are \code{"parzen"} and \code{"trapezoid"}.
   	}
  \item{ntoorderminus}{a value of \code{0,1}, or \code{2} such that the centering constant will retain terms of order \code{n^(-ntoorderminus)}.  Enter \code{0} for the moderate-p GCT, and enter \code{2} for the large-p GCT.  A value of \code{1} may be entered to retain only terms which are O(1/n), appropriate for a size of \code{p} between moderate and large.
  	} 
}

\value{
\item{TSvalue}{the unstudentized test statistic.
	}
\item{center}{the centering constant for studentizing the test statistic.
	}
\item{var}{the estimated variance of the unstudentized test statistic.
	}
\item{T}{the studentized test statistic.
	}
\item{smoother}{the choice of smoother used.
	}
\item{pvalue}{the p-value.
	}
}

\references{

	Gregory, K., Carroll, R. J., Baladandayuthapani, V. and Lahiri, S. (2015). A two-sample test for equality of means in high dimension. 
		\emph{Journal of the American Statistician}, to appear
}
\author{
Karl Gregory   \email{kgregory@mail.uni-mannheim.de}, \url{http://www.stat.tamu.edu/~kbgregory}.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{GCT.test.missing}}
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
	
	
	GCT.test(X,Y)
	
	}	
}

\keyword{htest }
\keyword{multivariate }% __ONLY ONE__ keyword per line
