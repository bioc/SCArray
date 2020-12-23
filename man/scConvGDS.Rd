\name{scConvGDS}
\alias{scConvGDS}
\title{Create a GDS file}
\description{
    Creates a single-cell GDS file from an R object.
}
\usage{
scConvGDS(obj, outfn, save.sp=TRUE, type=c("float32", "float64", "int32"),
    compress="LZMA_RA", clean=TRUE, verbose=TRUE)
}
\arguments{
    \item{obj}{a dense/sparse matrix, DelayedMatrix, SummarizedExperiment or
        SingleCellExperiment}
    \item{outfn}{the output file name in GDS format}
    \item{save.sp}{if \code{TRUE}, save it to a sparse matrix in GDS; otherwise,
        store dense matrix}
    \item{type}{numeric data type in the output file}
    \item{compress}{the compression method, see \code{\link{add.gdsn}}}
    \item{clean}{TRUE}
    \item{verbose}{if \code{TRUE}, show information}
}
\value{
    Return the path of the output file.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{scOpen}}, \code{\link{scClose}}
}

\examples{
# load a SingleCellExperiment object
fn <- system.file("extdata", "LaMannoBrainSub.rds", package="SCArray")
sce <- readRDS(fn)

scConvGDS(sce, "test.gds")

# remove the temporary output file
unlink("test.gds")
}

\keyword{GDS}
\keyword{SingleCell}
