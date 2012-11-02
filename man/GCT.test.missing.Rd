\name{GCT.test.missing}
\alias{GCT.test.missing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generalized component test for missing data
}
\description{Performs the generalized component test from Gregory et al. (2014) for the equality of two \code{p} by \code{1} population mean vectors given samples of sizes \code{n} and \code{m} when there are missing data.
}

\usage{
GCT.test.missing(X, Y, r, smoother = "parzen", ntoorderminus = 2)
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
\item{pvalue}{the pvalue.
	}
\item{overallpctmiss}{ the overall proportion of values that are missing.
	}
\item{pctmissperX}{a vector of length \code{p} containing the proportion of missing values per component in sample one.
	}
\item{pctmissperY}{a vector of length \code{p} containing the proportion of missing values per component in sample two.
	}
	
}

\references{
	Gregory, K., Carroll, R. J., Baladandayuthapani, V. and Lahiri, S. (2015). A two-sample test for equality of means in high dimension. 
		\emph{Journal of the American Statistician}, to appear
}
\author{
Karl Gregory  \email{kgregory@mail.uni-mannheim.de}, \url{http://www.stat.tamu.edu/~kbgregory}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{GCT.test}}
}
\examples{

	data(chr1qseg)

	X <- chr1qseg$X
	Y <- chr1qseg$Y
	
	GCT.test.missing(X,Y,r=20,smoother="parzen")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{multivariate}
