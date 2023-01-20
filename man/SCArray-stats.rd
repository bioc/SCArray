\name{SCArray-stats}
\alias{SCArray-stats}

\alias{rowSums}
\alias{rowSums,SC_GDSMatrix-method}
\alias{colSums}
\alias{colSums,SC_GDSMatrix-method}
\alias{rowSums2}
\alias{rowSums2,SC_GDSMatrix-method}
\alias{colSums2}
\alias{colSums2,SC_GDSMatrix-method}

\alias{rowProds}
\alias{rowProds,SC_GDSMatrix-method}
\alias{colProds}
\alias{colProds,SC_GDSMatrix-method}

\alias{rowMeans}
\alias{rowMeans,SC_GDSMatrix-method}
\alias{colMeans}
\alias{colMeans,SC_GDSMatrix-method}
\alias{rowMeans2}
\alias{rowMeans2,SC_GDSMatrix-method}
\alias{colMeans2}
\alias{colMeans2,SC_GDSMatrix-method}
\alias{rowWeightedMeans}
\alias{rowWeightedMeans,SC_GDSMatrix-method}
\alias{colWeightedMeans}
\alias{colWeightedMeans,SC_GDSMatrix-method}

\alias{rowVars}
\alias{rowVars,SC_GDSMatrix-method}
\alias{colVars}
\alias{colVars,SC_GDSMatrix-method}
\alias{rowWeightedVars}
\alias{rowWeightedVars,SC_GDSMatrix-method}
\alias{colWeightedVars}
\alias{colWeightedVars,SC_GDSMatrix-method}

\alias{rowSds}
\alias{rowSds,SC_GDSMatrix-method}
\alias{colSds}
\alias{colSds,SC_GDSMatrix-method}
\alias{rowWeightedSds}
\alias{rowWeightedSds,SC_GDSMatrix-method}
\alias{colWeightedSds}
\alias{colWeightedSds,SC_GDSMatrix-method}

\alias{rowMins}
\alias{rowMins,SC_GDSMatrix-method}
\alias{colMins}
\alias{colMins,SC_GDSMatrix-method}
\alias{rowMaxs}
\alias{rowMaxs,SC_GDSMatrix-method}
\alias{colMaxs}
\alias{colMaxs,SC_GDSMatrix-method}
\alias{rowRanges}
\alias{rowRanges,SC_GDSMatrix-method}
\alias{colRanges}
\alias{colRanges,SC_GDSMatrix-method}


\title{SC_GDSMatrix row/column summarization}
\description{
    The row/column summarization methods for the SC_GDSMatrix matrix,
extending the S4 methods in the \pkg{DelayedArray} and \pkg{DelayedMatrixStats}
packages.
}

\usage{
\S4method{rowSums}{SC_GDSMatrix}(x, na.rm=FALSE, dims=1)
\S4method{colSums}{SC_GDSMatrix}(x, na.rm=FALSE, dims=1)
\S4method{rowSums2}{SC_GDSMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
\S4method{colSums2}{SC_GDSMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)

\S4method{rowProds}{SC_GDSMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
\S4method{colProds}{SC_GDSMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)

\S4method{rowMeans}{SC_GDSMatrix}(x, na.rm=FALSE, dims=1)
\S4method{colMeans}{SC_GDSMatrix}(x, na.rm=FALSE, dims=1)
\S4method{rowMeans2}{SC_GDSMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
\S4method{colMeans2}{SC_GDSMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
\S4method{rowWeightedMeans}{SC_GDSMatrix}(x, w=NULL, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
\S4method{colWeightedMeans}{SC_GDSMatrix}(x, w=NULL, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)

\S4method{rowVars}{SC_GDSMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL, ..., useNames=NA)
\S4method{colVars}{SC_GDSMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL, ..., useNames=NA)
\S4method{rowWeightedVars}{SC_GDSMatrix}(x, w=NULL, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
\S4method{colWeightedVars}{SC_GDSMatrix}(x, w=NULL, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)

\S4method{rowSds}{SC_GDSMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL, ..., useNames=NA)
\S4method{colSds}{SC_GDSMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL, ..., useNames=NA)
\S4method{rowWeightedSds}{SC_GDSMatrix}(x, w=NULL, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
\S4method{colWeightedSds}{SC_GDSMatrix}(x, w=NULL, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)

\S4method{rowMins}{SC_GDSMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE)
\S4method{colMins}{SC_GDSMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE)
\S4method{rowMaxs}{SC_GDSMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE)
\S4method{colMaxs}{SC_GDSMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE)
\S4method{rowRanges}{SC_GDSMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE)
\S4method{colRanges}{SC_GDSMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE)
}

\arguments{
    \item{x}{A \link{SC_GDSMatrix} object (inherited from DelayedMatrix)}
    \item{dims}{not used, it should be 1}
    \item{rows, cols}{specify the subset of rows (and/or columns) to operate
        over; if \code{NULL}, no subsetting}
    \item{na.rm}{if \code{TRUE}, missing values (NaN and NA) will be removed}
    \item{w}{\code{NULL} or a numeric vector for weights}
    \item{center}{\code{NULL}, or a vector of pre-calculated row (column) means}
    \item{useNames}{if \code{TRUE}, the name attributes of result are set}
    \item{...}{additional arguments passed to specific methods}
}

\details{
    All these operations are block-processed according to the data stored in
the GDS file.
}

\seealso{
    \itemize{
        \item The \pkg{DelayedMatrixStats} package for more row/column
            summarization methods for \link{DelayedMatrix} objects.
        \item \link{DelayedMatrix-utils} for other common operations on
            \link{DelayedMatrix} objects.
        \item \link{DelayedMatrix} objects.
        \item \link[base]{matrix} objects in base R.
    }
}

\examples{
1
}

\keyword{methods}
\keyword{GDS}
