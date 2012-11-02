\name{build2popData}
\alias{build2popData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Data simulator for the high-dimensional two-sample setting
}
\description{ Generates data from two populations with user-specified mean vectors, covariance structure, sample sizes, and dimension of each observation.
}
\usage{
build2popData(n, m, p, muX, muY, dep, commoncov, VarScaleY, 
 ARMAparms, LRparm,S = 1, innov = function(n, ...) rnorm(n, 0, 1),
 heteroscedastic = FALSE, het.diag)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{number of observations in sample one.
}
  \item{m}{number of observations in sample two.
}
  \item{p}{number of components in each observation.
}
  \item{muX}{\code{p} by \code{1} vector of component means for population one.
}
  \item{muY}{\code{p} by \code{1} vector of component means for population two.
}
  \item{dep}{dependence structure among the \code{p} components for both populations. Choices are \code{'IND'}, \code{'ARMA'}, and \code{'LR'}.
}
  \item{commoncov}{a logical indicating whether populations one and two will have equal covariance matrices.  If \code{FALSE}, then the innovations used in generating the population two data will be scaled by the square root of the value specifed in \code{VarScaleY}.
}
  \item{VarScaleY}{constant by which innovations are scaled in generating observations for sample two when \code{commoncov=FALSE}.
}
  \item{ARMAparms}{a list of the form \code{list(coefs=list(ma= , ar=))} specifying the AR and MA parameters to be used when \code{dep='ARMA'}.
}
  \item{LRparm}{value of the LR dependence parameter to be used when \code{dep='LR'}.
}
  \item{S}{the number of data sets to simulate.
}
  \item{innov}{a function used to generate the innovations, such as \code{innov=function(n,...) rnorm(n,0,1)}.
}
  \item{heteroscedastic}{a logical indicating whether the components will be scaled by the entries in the diagonal matrix specified by \code{het.diag=}.
}
  \item{het.diag}{a \code{p} by \code{p} diagonal matrix, where the entries on the diagonal will be used to scale the component standard deviations.
}
}
\value{
A list of length \code{S} of lists, each containing
\item{X}{the \code{n} by \code{p} data matrix for sample one.}
\item{Y}{the \code{m} by \code{p} data matrix for sample two.}
\item{n}{the number of observations in sample one.}
\item{m}{the number of observation in sample two.}
\item{p}{the number of components in each observation.}
\item{muX}{the mean vector for population one.}
\item{muY}{the mean vector for population two.}
\item{dep}{the dependence structure chosen for \code{dep}.}
\item{commoncov}{the value of \code{commoncov} used in the function call.}
\item{VarScaleY}{the scalar by which the variance of the population two data is scaled.}
\item{ARMAparms}{the list containing the specified ARMA parameters.}
\item{LRparm}{the long-range dependence parameter.}
\item{S}{the number of simulated data sets.}
\item{innov}{the function chosen to generate the innovations.}
\item{heteroscedastic}{logical indicating whether \code{heteroscedastic=TRUE} in the function call.}
\item{het.diag}{the value of \code{het.diag} in the function call.}
}

\references{Hall, P. Jing, B. Y. and Lahiri, S. N. (1998). On the sampling window method for long-range dependent data.
			\emph{Statistica Sinica} \bold{8},1189--1204
}
\author{
Karl Gregory  \email{kgregory@mail.uni-mannheim.de}, \url{http://www.stat.tamu.edu/~kbgregory}.
}

\examples{
\dontrun{
DATA <-build2popData(	
	n = 15,
	m = 20,
	p = 500,
	muX = rep(0,500),
	muY = rep(0,500),
	commoncov = FALSE,
	VarscaleY = 2,
	dep = "ARMA",
	ARMAparms = list(coefs=list(ma=c(.2,.3) , ar=c(.4,-.1))),
	LRparm = .75,
	S = 25,
	innov = function(n,...) rnorm(n,0,1),
	heteroscedastic=TRUE,
	het.diag = diag(.1 + rexp(500,1/2))
	)
}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
\keyword{multivariate}% __ONLY ONE__ keyword per line
