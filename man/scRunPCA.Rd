\name{scRunPCA}
\alias{scRunPCA}
\alias{runPCA}
\alias{runPCA,SC_GDSMatrix-method}
\title{Perform PCA on SC_GDSMatrix and expression data}
\description{
    Perform a Principal Components Analysis (PCA) on cells in the
SingleCellExperiment object.
}
\usage{
scRunPCA(sce, ncomponents=50, ntop=500, subset_row=NULL, scale=FALSE,
    altexp=NULL, name="PCA", exprs_values="logcounts", dimred=NULL,
    n_dimred=NULL, BSPARAM=NULL, BPPARAM=SerialParam(), verbose=TRUE)

\S4method{runPCA}{SC_GDSMatrix}(x, rank, center=TRUE, scale=FALSE, get.rotation=TRUE,
    get.pcs=TRUE, ...)
}
\arguments{
    \item{sce}{a SingleCellExperiment or SummarizedExperiment object}
    \item{x}{a SC_GDSMatrix object}
    \item{ncomponents, rank}{# of calculated principal components}
    \item{ntop}{# of features with the highest variances to use for PCA}
    \item{subset_row}{specifying the subset of features to use}
    \item{center}{if \code{TRUE}, expression values will be be centered}
    \item{scale}{if \code{TRUE}, expression values will be be standardized}
    \item{altexp}{String or integer scalar specifying an alternative experiment
        containing the input data}
    \item{name}{the name to be used to store the result in \code{reducedDims}}
    \item{exprs_values}{the assay name containing the expression values}
    \item{dimred}{String or integer scalar specifying the existing
        dimensionality reduction results to use}
    \item{n_dimred}{Integer scalar or vector specifying the dimensions to use
        if \code{dimred} is specified}
    \item{BSPARAM}{A BiocSingularParam object specifying which algorithm to be
        used in \code{runPCA} in the BiocSingular package}
    \item{BPPARAM}{A BiocParallelParam object for parallelized calculation}
    \item{get.rotation}{if \code{TRUE}, return rotation vectors}
    \item{get.pcs}{if \code{TRUE}, return principal component scores}
    \item{verbose}{if TRUE, show information}
    \item{...}{For \code{runPCA}, this contains further arguments to pass
          to \code{runSVD}, including \code{BSPARAM} to specify the algorithm
          that should be used, and \code{BPPARAM} to control parallelization.}
}
\details{
    The function \code{runPCA()} simply calls \code{\link{runSVD}} and converts
the results into a format similar to that returned by \code{\link{prcomp}}.

\code{BSPARAM} can be one of
\describe{
    \item{\code{ExactParam()}:}{exact SVD with \code{\link{runExactSVD}}.}
    \item{\code{IrlbaParam()}:}{approximate SVD with \pkg{irlba} via
        \code{\link{runIrlbaSVD}}.}
    \item{\code{RandomParam()}:}{approximate SVD with \pkg{rsvd} via
        \code{\link{runRandomSVD}}.}
    \item{\code{FastAutoParam()}:}{fast approximate SVD, chosen based on
        the matrix representation.}
}
\code{fold=1} in BiocSingularParam is used for the situation that the
covariance matrix is relatively small, and running SVD on the small covariance
matrix can be more effecient. When \code{fold=Inf}, running SVD on the matrix
directly and will read the matrix multiple times. If it is a file-based matrix,
\code{fold=Inf} could be slow.
}
\value{
    Returns a SingleCellExperiment object containing the PC coordinate matrix
in \code{reducedDims(..., name)}. The attributes of the PC coordinate matrix
have "percentVar", "varExplained" and "rotation" (see \code{scater::runPCA}
for more details).
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{runSVD}} for the underlying SVD function.

    \code{?\linkS4class{BiocSingularParam}} for the SVD algorithm choices.

    \code{\link[BiocSingular]{runPCA}}.
}

\examples{
library(BiocSingular)

# a GDS file for SingleCellExperiment
fn <- system.file("extdata", "example.gds", package="SCArray")

x <- scArray(fn, "counts")
x <- x[1:200, ]
x

pc <- runPCA(x, BSPARAM=ExactParam(fold=1))  # using covariance matrix
str(pc)

rm(x)
}

\keyword{PCA}
\keyword{GDS}
\keyword{SingleCell}
