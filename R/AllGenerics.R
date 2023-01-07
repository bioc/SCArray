#######################################################################
#
# Package name: SCArray
#
# Description:
#     Large-scale single-cell RNA-seq data manipulation with GDS files
#
# Copyright (C) 2021-2023    Xiuwen Zheng (@AbbVie-ComputationalGenomics)
# License: GPL-3
#


#######################################################################

# For internal use only
setMethod("matrixClass", "SC_GDSArray", function(x) "SC_GDSMatrix")

# Automatic coercion method from GDSArray to GDSMatrix (muted for
# higher dimensions) this function works only when GDSArray is
# 2-dimensional, otherwise it fails.

setAs("SC_GDSArray", "SC_GDSMatrix", function(from) new("SC_GDSMatrix", from))
setAs("SC_GDSMatrix", "SC_GDSArray", function(from) from)
setAs("ANY", "SC_GDSMatrix",
    function(from) as(as(from, "SC_GDSArray"), "SC_GDSMatrix"))

# setMethod(
#     "DelayedArray", "SCArraySeed",
#     function(seed) new_DelayedArray(seed, Class="SC_GDSArray")
# )



#######################################################################

x_msg <- function(msg)
{
    if (getOption("scarray_debug", FALSE))
        message(msg)
}

x_check <- function(x, msg)
{
    if (getOption("SCArray_Debug", FALSE))
        message(msg)
    DelayedArray:::.get_ans_type(x, must.be.numeric=TRUE)
    stopifnot(length(dim(x)) == 2L)
    invisible()
}

x_subset <- function(x, rows, cols)
{
    if (!is.null(rows) && !is.null(cols))
    {
        x <- x[rows, cols, drop = FALSE]
    } else if (!is.null(rows))
    {
        x <- x[rows, , drop = FALSE]
    } else if (!is.null(cols))
    {
        x <- x[, cols, drop = FALSE]
    }
    x
}

x_nperm <- function(x)
{
    ans <- if (is(x, "DelayedAperm")) 1L else 0L
    if (is(x, "DelayedUnaryOp"))
        return(ans + x_nperm(x@seed))
    if (is(x, "DelayedNaryOp"))
        x <- x@seeds
    if (is.list(x) && !is.array(x))
    {
        for (y in x) ans <- ans + x_nperm(y)
    }
    ans
}

x_type <- function(x)
{
    if (nseed(x)==1L && is(seed(x), "SCArraySeed"))
    {
        # whether transposed or not
        if (x_nperm(x) %% 2L == 0L) 1L else 2L
    } else
        3L
}


################

.x_row_sums <- function(x, na.rm, ...)
{
    # block read
    blockReduce(function(bk, v, na.rm) {
        .Call(c_rowSums_update, bk, v, na.rm)
    }, x, init=double(nrow(x)), grid=colAutoGrid(x), na.rm=na.rm, ...)
}

x_rowSums <- function(x, na.rm=FALSE, dims=1)
{
    x_check(x, "Calling SCArray:::x_rowSums() ...")
    stopifnot(identical(dims, 1))
    if (x_type(x) == 1L)
    {
        # output
        .x_row_sums(x, na.rm)
    } else {
        x_msg("Calling DelayedArray::rowSums() ...")
        callNextMethod()
    }
}

x_colSums <- function(x, na.rm=FALSE, dims=1)
{
    x_check(x, "Calling SCArray:::x_colSums() ...")
    stopifnot(identical(dims, 1))
    if (x_type(x) != 2L)
    {
        x_msg("Calling DelayedArray::colSums() ...")
        callNextMethod()
    } else {
        x_msg("Calling DelayedArray::.x_row_sums() ...")
        .x_row_sums(t(x), na.rm)
    }
}

x_rowSums2 <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
{
    x_check(x, "Calling SCArray:::x_rowSums2() ...")
    if (x_type(x) == 1L)
    {
        # output
        x <- x_subset(x, rows, cols)
        v <- .x_row_sums(x, na.rm, ...)
        if (isTRUE(useNames)) names(v) <- rownames(x)
        v
    } else {
        x_msg("Calling DelayedMatrixStats::rowSums2() ...")
        callNextMethod()
    }
}

x_colSums2 <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
{
    x_check(x, "Calling SCArray:::x_colSums2() ...")
    if (x_type(x) != 2L)
    {
        x_msg("Calling DelayedMatrixStats::colSums2() ...")
        callNextMethod()
    } else {
        x_msg("Calling DelayedArray::.x_row_sums() ...")
        x <- x_subset(x, rows, cols)
        v <- .x_row_sums(t(x), na.rm, ...)
        if (isTRUE(useNames)) names(v) <- colnames(x)
        v
    }
}

setMethod("rowSums", "SC_GDSMatrix", x_rowSums)
setMethod("colSums", "SC_GDSMatrix", x_colSums)
setMethod("rowSums2", "SC_GDSMatrix", x_rowSums2)
setMethod("colSums2", "SC_GDSMatrix", x_colSums2)


################

.x_row_prods <- function(x, na.rm, ...)
{
    # block read
    blockReduce(function(bk, v, na.rm) {
        .Call(c_rowProds_update, bk, v, na.rm)
    }, x, init=rep(1, nrow(x)), grid=colAutoGrid(x), na.rm=na.rm, ...)
}

x_rowProds <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
{
    x_check(x, "Calling SCArray:::x_rowProds() ...")
    if (x_type(x) == 1L)
    {
        # output
        x <- x_subset(x, rows, cols)
        v <- .x_row_prods(x, na.rm, ...)
        if (isTRUE(useNames)) names(v) <- rownames(x)
        v
    } else {
        x_msg("Calling DelayedMatrixStats::rowProds() ...")
        callNextMethod()
    }
}

x_colProds <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
{
    x_check(x, "Calling SCArray:::x_colProds() ...")
    if (x_type(x) != 2L)
    {
        x_msg("Calling DelayedMatrixStats::colProds() ...")
        callNextMethod()
    } else {
        x_msg("Calling DelayedArray::.x_row_prods() ...")
        x <- x_subset(x, rows, cols)
        v <- .x_row_prods(t(x), na.rm, ...)
        if (isTRUE(useNames)) names(v) <- colnames(x)
        v
    }
}

setMethod("rowProds", "SC_GDSMatrix", x_rowProds)
setMethod("colProds", "SC_GDSMatrix", x_colProds)


################

.x_row_means <- function(x, na.rm, ...)
{
    # block read
    rv <- blockReduce(function(bk, v, na.rm) {
        .Call(c_rowMeans_update, bk, v, na.rm)
    }, x, init=double(nrow(x)*2L), grid=colAutoGrid(x), na.rm=na.rm, ...)
    # finally
    .Call(c_rowMeans_final, rv)
}

x_rowMeans <- function(x, na.rm=FALSE, dims=1)
{
    x_check(x, "Calling SCArray:::x_rowMeans() ...")
    stopifnot(identical(dims, 1))
    if (x_type(x) == 1L)
    {
        # output
        .x_row_means(x, na.rm)
    } else {
        x_msg("Calling DelayedArray::rowMeans() ...")
        callNextMethod()
    }
}

x_colMeans <- function(x, na.rm=FALSE, dims=1)
{
    x_check(x, "Calling SCArray:::x_colMeans() ...")
    stopifnot(identical(dims, 1))
    if (x_type(x) != 2L)
    {
        x_msg("Calling DelayedArray::colMeans() ...")
        callNextMethod()
    } else {
        x_msg("Calling DelayedArray::.x_row_means() ...")
        .x_row_means(t(x), na.rm)
    }
}

x_rowMeans2 <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
{
    x_check(x, "Calling SCArray:::x_rowMeans2() ...")
    if (x_type(x) == 1L)
    {
        # output
        x <- x_subset(x, rows, cols)
        v <- .x_row_means(x, na.rm, ...)
        if (isTRUE(useNames)) names(v) <- rownames(x)
        v
    } else {
        x_msg("Calling DelayedMatrixStats::rowMeans2() ...")
        callNextMethod()
    }
}

x_colMeans2 <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
{
    x_check(x, "Calling SCArray:::x_colMeans2() ...")
    if (x_type(x) != 2L)
    {
        x_msg("Calling DelayedMatrixStats::colMeans2() ...")
        callNextMethod()
    } else {
        x_msg("Calling DelayedArray::.x_row_means() ...")
        x <- x_subset(x, rows, cols)
        v <- .x_row_means(t(x), na.rm, ...)
        if (isTRUE(useNames)) names(v) <- colnames(x)
        v
    }
}

setMethod("rowMeans", "SC_GDSMatrix", x_rowMeans)
setMethod("colMeans", "SC_GDSMatrix", x_colMeans)
setMethod("rowMeans2", "SC_GDSMatrix", x_rowMeans2)
setMethod("colMeans2", "SC_GDSMatrix", x_colMeans2)


################

.x_row_vars <- function(x, na.rm, center, ...)
{
    # check
    if (!is.null(center))
    {
        stopifnot(nrow(x) == length(center))
        if (is.integer(center)) center <- as.double(center)
    }
    # block read
    v <- blockReduce(function(bk, v, na.rm, center) {
        .Call(c_rowVars_update, bk, v, na.rm, center)
    }, x, init=double(nrow(x)*3L), grid=colAutoGrid(x),
        na.rm=na.rm, center=center, ...)
    # finally
    .Call(c_rowVars_final, v, center)
}

x_rowVars <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL, ...,
    useNames=NA)
{
    x_check(x, "Calling SCArray:::x_rowVars() ...")
    stopifnot(is.null(center) || is.numeric(center))
    if (x_type(x) == 1L)
    {
        x <- x_subset(x, rows, cols)
        v <- .x_row_vars(x, na.rm, center, ...)
        if (isTRUE(useNames)) names(v) <- rownames(x)
        v
    } else {
        x_msg("Calling DelayedMatrixStats::rowVars() ...")
        callNextMethod()
    }
}

x_colVars <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL, ...,
    useNames=NA)
{
    x_check(x, "Calling SCArray:::x_colVars() ...")
    if (x_type(x) != 2L)
    {
        x_msg("Calling DelayedMatrixStats::colVars() ...")
        callNextMethod()
    } else {
        x_msg("Calling DelayedArray::.x_row_vars() ...")
        x <- x_subset(x, rows, cols)
        v <- .x_row_vars(t(x), na.rm, center, ...)
        if (isTRUE(useNames)) names(v) <- colnames(x)
        v
    }
}

setMethod("rowVars", "SC_GDSMatrix", x_rowVars)
setMethod("colVars", "SC_GDSMatrix", x_colVars)


################

x_rowSds <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL, ...,
    useNames=NA)
{
    x_check(x, "Calling SCArray:::x_rowSds() ...")
    v <- x_rowVars(x, rows, cols, na.rm, center, ..., useNames=useNames)
    sqrt(v)
}

x_colSds <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL, ...,
    useNames=NA)
{
    x_check(x, "Calling SCArray:::x_colSds() ...")
    v <- x_colVars(x, rows, cols, na.rm, center, ..., useNames=useNames)
    sqrt(v)
}

setMethod("rowSds", "SC_GDSMatrix", x_rowSds)
setMethod("colSds", "SC_GDSMatrix", x_colSds)



################

# \S4method{rowMins}{DelayedMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE)
# \S4method{colMins}{DelayedMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE)

# \S4method{rowMaxs}{DelayedMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE)
# \S4method{colMaxs}{DelayedMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE)

# \S4method{rowRanges}{DelayedMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE)
# \S4method{colRanges}{DelayedMatrix}(x, rows=NULL, cols=NULL, na.rm=FALSE)




