\name{scOpen}
\alias{scOpen}
\title{Open a Single-cell GDS File}
\description{
    Opens a single-cell GDS file.
}
\usage{
scOpen(gdsfn, readonly=TRUE, allow.duplicate=FALSE)
}
\arguments{
    \item{gdsfn}{the input file name}
    \item{readonly}{whether read-only or not}
    \item{allow.duplicate}{if \code{TRUE}, it is allowed to open a GDS file
        with read-only mode when it has been opened in the same R session}
}
\value{
    Return an object of class \code{SCArrayFileClass} inherited from
\code{\link{gds.class}}.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{scClose}}
}

\examples{
# check
cat("OK\n")
}

\keyword{GDS}
\keyword{SingleCell}
