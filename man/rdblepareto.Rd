\name{rdblepareto}
\alias{rdblepareto}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The double pareto distribution
%%  ~~function to do ... ~~
}
\description{Random number generation from the double Pareto distribution with parameters \code{shape} and \code{scale}.}
\usage{
rdblepareto(n, shape, scale)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ number of observations to generate. 
}
  \item{shape}{shape parameter, must be strictly greater than zero.
}
  \item{scale}{scale parameter, must be strictly greater than zero.
}
}

\value{ a vector of length \code{n} of realizations from the double Pareto distribution with shape parameter \code{shape} and scale parameter \code{scale}.
}

\author{
Karl Gregory   \email{kgregory@mail.uni-mannheim.de}, \url{http://www.stat.tamu.edu/~kbgregory}.
}

\examples{

## The function is currently defined as
function (n, shape, scale) 
{
    u <- runif(n, 0, 1)
    x <- scale * (1 - u)^(-1/shape) - scale
    y <- sample(c(-1, 1), n, replace = TRUE) * x
    return(y)
  }


x <- rdblepareto(500,1.5,1)

plot(density(x))
}
\keyword{distribution}
\keyword{datagen}
