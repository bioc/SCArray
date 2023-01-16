\name{scRunPCA}
\alias{scRunPCA}
\title{Perform PCA on expression data}
\description{
    Perform a Principal Components Analysis (PCA) on cells in the
SingleCellExperiment object.
}
\usage{
scRunPCA(x, ncomponents=50, ntop=500, subset_row=NULL, scale=FALSE,
    altexp=NULL, name="PCA", exprs_values="logcounts", dimred=NULL,
    n_dimred=NULL, BSPARAM=NULL, BPPARAM=SerialParam(), verbose=TRUE)
}
\arguments{
    \item{x}{a SingleCellExperiment or SummarizedExperiment object}
    \item{ncomponents}{# of calculated principal components}
    \item{ntop}{# of features with the highest variances to use for PCA}
    \item{subset_row}{specifying the subset of features to use}
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
    \item{verbose}{if TRUE, show information}
}
\value{
    Return an object with the class \code{SC_GDSMatrix} or \code{SC_GDSArray}.
}

\author{Xiuwen Zheng}
\seealso{
    \code{\link{scExperiment}}
}

\examples{
1
}

\keyword{PCA}
\keyword{GDS}
\keyword{SingleCell}
