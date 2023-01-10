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
    if (getOption("scarray_debug", FALSE))
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

# Return 1 for SCArraySeed, 2L for transposed SCArraySeed, 3L for others
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
    blockReduce(function(bk, v, na.rm) {
        .Call(c_rowSums, bk, v, na.rm)
    }, x, double(nrow(x)), grid=colAutoGrid(x), as.sparse=NA, na.rm=na.rm, ...)
}

x_rowSums <- function(x, na.rm=FALSE, dims=1)
{
    x_check(x, "Calling SCArray:::x_rowSums() ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
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
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
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
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    if (x_type(x) == 1L)
    {
        x <- x_subset(x, rows, cols)
        v <- .x_row_sums(x, na.rm, ...)
        if (isTRUE(useNames)) names(v) <- rownames(x)
        v  # output
    } else {
        x_msg("Calling DelayedMatrixStats::rowSums2() ...")
        callNextMethod()
    }
}

x_colSums2 <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
{
    x_check(x, "Calling SCArray:::x_colSums2() ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    if (x_type(x) != 2L)
    {
        x_msg("Calling DelayedMatrixStats::colSums2() ...")
        callNextMethod()
    } else {
        x_msg("Calling DelayedArray::.x_row_sums() ...")
        x <- x_subset(x, rows, cols)
        v <- .x_row_sums(t(x), na.rm, ...)
        if (isTRUE(useNames)) names(v) <- colnames(x)
        v  # output
    }
}

setMethod("rowSums", "SC_GDSMatrix", x_rowSums)
setMethod("colSums", "SC_GDSMatrix", x_colSums)
setMethod("rowSums2", "SC_GDSMatrix", x_rowSums2)
setMethod("colSums2", "SC_GDSMatrix", x_colSums2)


################

.x_row_prods <- function(x, na.rm, ...)
{
    blockReduce(function(bk, v, na.rm) {
        .Call(c_rowProds, bk, v, na.rm)
    }, x, init=rep(1, nrow(x)), grid=colAutoGrid(x), na.rm=na.rm, ...)
}

x_rowProds <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
{
    x_check(x, "Calling SCArray:::x_rowProds() ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    if (x_type(x) == 1L)
    {
        x <- x_subset(x, rows, cols)
        v <- .x_row_prods(x, na.rm, ...)
        if (isTRUE(useNames)) names(v) <- rownames(x)
        v  # output
    } else {
        x_msg("Calling DelayedMatrixStats::rowProds() ...")
        callNextMethod()
    }
}

x_colProds <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
{
    x_check(x, "Calling SCArray:::x_colProds() ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    if (x_type(x) != 2L)
    {
        x_msg("Calling DelayedMatrixStats::colProds() ...")
        callNextMethod()
    } else {
        x_msg("Calling DelayedArray::.x_row_prods() ...")
        x <- x_subset(x, rows, cols)
        v <- .x_row_prods(t(x), na.rm, ...)
        if (isTRUE(useNames)) names(v) <- colnames(x)
        v  # output
    }
}

setMethod("rowProds", "SC_GDSMatrix", x_rowProds)
setMethod("colProds", "SC_GDSMatrix", x_colProds)


################

.x_row_means <- function(x, na.rm, ...)
{
    # block read
    rv <- blockReduce(function(bk, v, na.rm) {
        .Call(c_rowMeans, bk, v, na.rm)
    }, x, double(nrow(x)*2L), grid=colAutoGrid(x), as.sparse=NA, na.rm=na.rm, ...)
    # finally
    .Call(c_rowMeans_final, rv)
}

.x_col_means <- function(x, na.rm, ...)
{
    unlist(blockApply(x, function(bk, na.rm) {
        .Call(c_colMeans, bk, na.rm)
    }, grid=colAutoGrid(x), as.sparse=NA, na.rm=na.rm, ...))
}

x_rowMeans <- function(x, na.rm=FALSE, dims=1)
{
    x_check(x, "Calling SCArray:::x_rowMeans() ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    stopifnot(identical(dims, 1))
    # output
    switch(x_type(x),
        .x_row_means(x, na.rm),     # 1
        .x_col_means(t(x), na.rm),  # 2
        callNextMethod()            # 3
    )
}

x_colMeans <- function(x, na.rm=FALSE, dims=1)
{
    x_check(x, "Calling SCArray:::x_colMeans() ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    stopifnot(identical(dims, 1))
    switch(x_type(x),
        .x_col_means(x, na.rm),     # 1
        .x_row_means(t(x), na.rm),  # 2
        callNextMethod()            # 3
    )
}

x_rowMeans2 <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, ..., useNames=NA)
{
    x_check(x, "Calling SCArray:::x_rowMeans2() ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    k <- x_type(x)
    if (k < 3L)
    {
        x <- x_subset(x, rows, cols)
        if (k == 1L)
            v <- .x_row_means(x, na.rm, ...)
        else
            v <- .x_col_means(t(x), na.rm, ...)
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
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    k <- x_type(x)
    if (k < 3L)
    {
        x <- x_subset(x, rows, cols)
        if (k == 1L)
            v <- .x_col_means(x, na.rm, ...)
        else
            v <- .x_row_means(t(x), na.rm, ...)
        if (isTRUE(useNames)) names(v) <- colnames(x)
        v
    } else {
        x_msg("Calling DelayedMatrixStats::colMeans2() ...")
        callNextMethod()
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
        .Call(c_rowVars, bk, v, na.rm, center)
    }, x, init=double(nrow(x)*3L), grid=colAutoGrid(x), as.sparse=NA,
        na.rm=na.rm, center=center, ...)
    # finally
    .Call(c_rowVars_final, v, center)
}

.x_col_vars <- function(x, na.rm, center, ...)
{
    # check
    if (!is.null(center))
    {
        stopifnot(ncol(x) == length(center))
        if (is.integer(center)) center <- as.double(center)
    }
    # block read
    unlist(blockApply(x, function(bk, na.rm, center) {
        .Call(c_colVars, bk, na.rm, center)
    }, grid=colAutoGrid(x), as.sparse=NA, na.rm=na.rm, center=center, ...))
}

x_rowVars <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, center=NULL, ...,
    useNames=NA)
{
    x_check(x, "Calling SCArray:::x_rowVars() ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    stopifnot(is.null(center) || is.numeric(center))
    k <- x_type(x)
    if (k < 3L)
    {
        x <- x_subset(x, rows, cols)
        if (k == 1L)
            v <- .x_row_vars(x, na.rm, center, ...)
        else
            v <- .x_col_vars(t(x), na.rm, center, ...)
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
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    stopifnot(is.null(center) || is.numeric(center))
    k <- x_type(x)
    if (k < 3L)
    {
        x <- x_subset(x, rows, cols)
        if (k == 1L)
            v <- .x_col_vars(x, na.rm, center, ...)
        else
            v <- .x_row_vars(t(x), na.rm, center, ...)
        if (isTRUE(useNames)) names(v) <- colnames(x)
        v
    } else {
        x_msg("Calling DelayedMatrixStats::colVars() ...")
        callNextMethod()
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

.x_row_mins <- function(x, na.rm, ...)
{
    blockReduce(function(bk, v, na.rm) {
        .Call(c_rowMins, bk, v, na.rm)
    }, x, rep(Inf, nrow(x)), grid=colAutoGrid(x), as.sparse=NA, na.rm=na.rm, ...)
}

.x_col_mins <- function(x, na.rm, ...)
{
    unlist(blockApply(x, function(bk, na.rm) {
        .Call(c_colMins, bk, na.rm)
    }, grid=colAutoGrid(x), as.sparse=NA, na.rm=na.rm, ...))
}

x_rowMins <- function(x, rows=NULL, cols=NULL, na.rm=FALSE)
{
    x_check(x, "Calling SCArray:::x_rowMins() ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    switch(x_type(x),
        .x_row_mins(x_subset(x, rows, cols), na.rm),     # 1
        .x_col_mins(t(x_subset(x, rows, cols)), na.rm),  # 2
        {   # 3
            x_msg("Calling DelayedArray::rowMins() ...")
            callNextMethod()
        })
}

x_colMins <- function(x, rows=NULL, cols=NULL, na.rm=FALSE)
{
    x_check(x, "Calling SCArray:::x_colMins() ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    switch(x_type(x),
        .x_col_mins(x_subset(x, rows, cols), na.rm),     # 1
        .x_row_mins(t(x_subset(x, rows, cols)), na.rm),  # 2
        {   # 3
            x_msg("Calling DelayedArray::colMins() ...")
            callNextMethod()
        })
}

setMethod("rowMins", "SC_GDSMatrix", x_rowMins)
setMethod("colMins", "SC_GDSMatrix", x_colMins)


################

.x_row_maxs <- function(x, na.rm, ...)
{
    blockReduce(function(bk, v, na.rm) {
        .Call(c_rowMaxs, bk, v, na.rm)
    }, x, rep(-Inf, nrow(x)), grid=colAutoGrid(x), as.sparse=NA, na.rm=na.rm, ...)
}

.x_col_maxs <- function(x, na.rm, ...)
{
    unlist(blockApply(x, function(bk, na.rm) {
        .Call(c_colMaxs, bk, na.rm)
    }, grid=colAutoGrid(x), as.sparse=NA, na.rm=na.rm, ...))
}

x_rowMaxs <- function(x, rows=NULL, cols=NULL, na.rm=FALSE)
{
    x_check(x, "Calling SCArray:::x_rowMaxs() ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    switch(x_type(x),
        .x_row_maxs(x_subset(x, rows, cols), na.rm),     # 1
        .x_col_maxs(t(x_subset(x, rows, cols)), na.rm),  # 2
        {   # 3
            x_msg("Calling DelayedArray::rowMaxs() ...")
            callNextMethod()
        })
}

x_colMaxs <- function(x, rows=NULL, cols=NULL, na.rm=FALSE)
{
    x_check(x, "Calling SCArray:::x_colMaxs() ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    switch(x_type(x),
        .x_col_maxs(x_subset(x, rows, cols), na.rm),     # 1
        .x_row_maxs(t(x_subset(x, rows, cols)), na.rm),  # 2
        {   # 3
            x_msg("Calling DelayedArray::colMaxs() ...")
            callNextMethod()
        })
}

setMethod("rowMaxs", "SC_GDSMatrix", x_rowMaxs)
setMethod("colMaxs", "SC_GDSMatrix", x_colMaxs)


################

.x_row_ranges <- function(x, na.rm, ...)
{
    blockReduce(function(bk, v, na.rm) {
        .Call(c_rowRanges, bk, v, na.rm)
    }, x, init=t(matrix(c(Inf, -Inf), nrow=2L, ncol=nrow(x))),
        grid=colAutoGrid(x), as.sparse=NA, na.rm=na.rm, ...)
}

.x_col_ranges <- function(x, na.rm, ...)
{
    do.call(rbind, blockApply(x, function(bk, na.rm) {
        .Call(c_colRanges, bk, na.rm)
    }, grid=colAutoGrid(x), as.sparse=NA, na.rm=na.rm, ...))
}

x_rowRanges <- function(x, rows=NULL, cols=NULL, na.rm=FALSE)
{
    x_check(x, "Calling SCArray:::x_rowRanges() ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    switch(x_type(x),
        .x_row_ranges(x_subset(x, rows, cols), na.rm),     # 1
        .x_col_ranges(t(x_subset(x, rows, cols)), na.rm),  # 2
        {   # 3
            x_msg("Calling DelayedArray::rowRanges() ...")
            callNextMethod()
        })
}

x_colRanges <- function(x, rows=NULL, cols=NULL, na.rm=FALSE)
{
    x_check(x, "Calling SCArray:::x_colRanges() ...")
    stopifnot(is.logical(na.rm), length(na.rm)==1L)
    switch(x_type(x),
        .x_col_ranges(x_subset(x, rows, cols), na.rm),     # 1
        .x_row_ranges(t(x_subset(x, rows, cols)), na.rm),  # 2
        {   # 3
            x_msg("Calling DelayedArray::colRanges() ...")
            callNextMethod()
        })
}

setMethod("rowRanges", "SC_GDSMatrix", x_rowRanges)
setMethod("colRanges", "SC_GDSMatrix", x_colRanges)





