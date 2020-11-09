\name{scOpen}
\alias{scOpen}
\title{Open a Single-cell GDS File}
\description{
    Opens a single-cell GDS file.
}
\usage{
scOpen(gdsfn, readonly=TRUE, allow.duplicate=TRUE)
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
    \code{\link{scClose}}, \code{\link{scArray}}
}

\examples{
# a GDS file for SingleCellExperiment
fn <- system.file("extdata", "LaMannoBrainData.gds", package="SCArray")

# open the GDS file
(f <- scOpen(fn))

# read a GDS file
cell.id <- read.gdsn(index.gdsn(f, "feature.id"))
samp.id <- read.gdsn(index.gdsn(f, "sample.id"))

# get a DelayedArray object
(cnt <- scArray(f, "counts"))

scClose(f)
}

\keyword{GDS}
\keyword{SingleCell}
