\name{scExperiment}
\alias{scExperiment}
\title{Get a SummarizedExperiment}
\description{
    Gets an instance of SingleCellExperiment or SummarizedExperiment.
}
\usage{
scExperiment(gdsfile, sc=TRUE, use.names=TRUE)
}
\arguments{
    \item{gdsfile}{a single-cell GDS object with class \code{SCArrayFileClass}}
    \item{sc}{if \code{TRUE}, return an instance of \code{SingleCellExperiment};
        otherwise an instance of \code{SummarizedExperiment}}
    \item{use.names}{if \code{TRUE}, load \code{dimnames} from 'feature.id' and
        'sample.id'}
}
\value{
    Return an instance of \code{\link{SingleCellExperiment}} or
\code{\link{SummarizedExperiment}}.
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
