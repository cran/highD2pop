\name{rgammashift}
\alias{rgammashift}
\title{The centered gamma distribution
}
\description{Random generation for the Gamma distribution with parameters \code{shape} and \code{scale} centered so as to have mean zero.}
\usage{
rgammashift(n, shape, scale)
}

\arguments{
  \item{n}{ number of observations to generate.
}
  \item{shape}{shape parameter, must be strictly greater than zero.
}
  \item{scale}{scale parameter, must be strictly greater than zero.
}
}

\value{ a vector of length \code{n} of realizations from the centered gamma distribution with shape parameter \code{shape} and scale parameter \code{scale}.

}
\author{
Karl Gregory  \email{kgregory@mail.uni-mannheim.de}, \url{http://www.stat.tamu.edu/~kbgregory}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{

## The function is currently defined as
function (n, shape, scale) 
{
    x <- rgamma(n = n, shape = shape, scale = scale) - shape*scale
    return(x)
}


x <- rgammashift(500,4,2)

plot(density(x))

}


\keyword{distribution}
\keyword{datagen}
