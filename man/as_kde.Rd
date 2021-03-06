% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distribution_functions.R
\name{as_kde}
\alias{as_kde}
\alias{as_kde.density}
\title{Coercion to class kde}
\usage{
as_kde(x, ...)

\method{as_kde}{density}(x, ...)
}
\arguments{
\item{x}{an object}

\item{...}{extra arguments for methods.  Currently not used.}
}
\value{
an object of class \code{"kde"}
}
\description{
A coercion function to class \code{"kde"}, currently the
only method is for objects of class \code{"density"} from
the \code{stats} package.
}
\examples{
set.seed(1234)
.x <- rnorm(5000)
dx <- density(.x)
oldPar <- par(mfrow = c(1,2), pty = "s", cex.main = 0.8,
  mar = c(3, 3, 1, 1) + 0.1, las = 1, mgp = c(1.5, 0.5, 0),
  cex.axis = 0.7, cex.lab = 0.8, tck = -0.015)
plot(dx)
plot(as_kde(dx))
par(oldPar)
rm(.x, oldPar)
}
