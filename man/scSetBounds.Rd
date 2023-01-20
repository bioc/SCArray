\name{scSetBounds}
\alias{scSetBounds}
\alias{scSetMax}
\alias{scSetMin}
\title{Set the bounds}
\description{
    Set the maximum and/or minimum on a GDS-specific DelayedArray.
}
\usage{
scSetMax(x, vmax)
scSetMin(x, vmin)
scSetBounds(x, vmin=NaN, vmax=NaN)
}
\arguments{
    \item{x}{a \code{SC_GDSArray} object}
    \item{vmax}{maximum}
    \item{vmin}{minimum}
}
\value{
    Return an object with the class \code{SC_GDSMatrix} or \code{SC_GDSArray}.
}

\author{Xiuwen Zheng}

\examples{
1
}

\keyword{PCA}
\keyword{GDS}
\keyword{SingleCell}
