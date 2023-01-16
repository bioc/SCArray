\name{scObj}
\alias{scObj}
\title{DelayedArray Object in GDS}
\description{
    Convert to SC_GDSArray/SC_GDSMatrix for utilizing seed-aware DelayedArray
functions.
}
\usage{
scObj(obj, verbose=FALSE)
}
\arguments{
    \item{obj}{a SummarizedExperiment, SingleCellExperiment or DelayedArray
        object}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    Return the object \code{obj} with the object class \code{DelayedArray}
replaced by the class \code{SC_GDSMatrix} or \code{SC_GDSArray}.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{scArray}}, \code{\link{scExperiment}}
}

\examples{
# a GDS file for SingleCellExperiment
fn <- system.file("extdata", "LaMannoBrainData.gds", package="SCArray")

cnt <- scArray(fn, "counts")
cnt  # it is SC_GDSMatrix

logcnt <- log2(cnt + 1)
logcnt  # it is a general DelayedMatrix

# convert to SC_GDSMatrix
lcnt <- scObj(logcnt)
lcnt
}

\keyword{GDS}
\keyword{SingleCell}
