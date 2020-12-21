\name{scClose}
\alias{scClose}
\title{Close the Single-cell GDS File}
\description{
    Closes a single-cell GDS file which is open.
}
\usage{
scClose(gdsfile)
}
\arguments{
    \item{gdsfile}{a single-cell GDS object with class \code{SCArrayFileClass}}
}
\value{
    None.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{scOpen}}
}

\examples{
# a GDS file for SingleCellExperiment
fn <- system.file("extdata", "LaMannoBrainData.gds", package="SCArray")

# open the GDS file
(f <- scOpen(fn))

# read
cell.id <- read.gdsn(index.gdsn(f, "feature.id"))
samp.id <- read.gdsn(index.gdsn(f, "sample.id"))

scClose(f)
}

\keyword{GDS}
\keyword{SingleCell}
