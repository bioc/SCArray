\name{scObj}
\alias{scObj}
\title{Seed-aware DelayedArray Object}
\description{
    Convert to SC_GDSArray for utilizing seed-aware DelayedArray functions.
}
\usage{
scObj(obj, verbose=TRUE)
}
\arguments{
    \item{obj}{a SummarizedExperiment, SingleCellExperiment or DelayedArray
        object}
    \item{verbose}{if TRUE, show information}
}
\value{
    Return an object with \code{SC_GDSMatrix}/\code{SC_GDSArray}.
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
