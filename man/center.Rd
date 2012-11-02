\name{center}
\alias{center}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The centering value for the test statistic
%%  ~~function to do ... ~~
}
\description{This is a helper function for the generalized component test from Gregory et al. (2014) which computes and returns the value \eqn{1 + \hat{c}_{jn}/n + \hat{d}_{jn}/n^2} for a component \eqn{j}. The \code{\link{GCT.test}} and \code{\link{GCT.sim}} average these values across \eqn{j=1,\dots,p} to obtain the centering \eqn{1 + \hat{a}_{n}/n + \hat{b}_{n}/n^2} for the generalized component test statistic.}
\usage{
center(xy,n,m,ntoorderminus=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xy}{a vector of length \code{n + m} containing the concatenation of the \code{n} realizations of a given component from the first sample and the \code{m} realizations of the same component from the second sample.
}
  \item{n}{the number of entries in \code{xy} belonging to the first sample.}

 \item{m}{the number of entries in \code{xy} belonging to the second sample.}

  \item{ntoorderminus}{the order of terms retained in \eqn{1 + \hat{c}_{jn}/n + \hat{d}_{jn}/n^2}.  If \code{ntoorderminus=0}, then \eqn{1} is returned. If \code{ntoorderminus=1} then \eqn{1 + \hat{c}_{jn}/n} is returned.  If \code{ntoorderminus=2} then \eqn{1 + \hat{c}_{jn}/n+\hat{d}_{jn}/n^2 } is returned.}

}

\value{A scalar equal to \eqn{1}, \eqn{1 + \hat{c}_{jn}/n} , or \eqn{1 + \hat{c}_{jn}/n+\hat{d}_{jn}/n^2 }, respectively, for \code{ntoorderminus=0,1,2}.
}
\references{
	Gregory, K., Carroll, R. J., Baladandayuthapani, V. and Lahiri, S. (2015). A two-sample test for equality of means in high dimension. 
		\emph{Journal of the American Statistician}, to appear}
\author{
Karl Gregory   \email{kgregory@mail.uni-mannheim.de}, \url{http://www.stat.tamu.edu/~kbgregory}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\keyword{teststatistic}
