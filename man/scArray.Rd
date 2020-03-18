\name{scArray}
\alias{scArray}
\title{Get an DelayedArray object}
\description{
    Gets an DelayedArray object from a single-cell GDS file.
}
\usage{
scArray(gdsfile, varname)
}
\arguments{
    \item{gdsfile}{a single-cell GDS object with class \code{SCArrayFileClass}}
    \item{varname}{character for the node name in the GDS file}
}
\value{
    Return an object of class \code{\link{DelayedArray}}.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{scOpen}}, \code{\link{scClose}}
}

\examples{
# check
cat("OK\n")
}

\keyword{GDS}
\keyword{SingleCell}
