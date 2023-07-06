#######################################################################
#
# Package name: SCArray
#
# Description:
#     Large-scale single-cell RNA-seq data manipulation with GDS files
#
# Copyright (C) 2023    Xiuwen Zheng (@AbbVie-ComputationalGenomics)
# License: GPL-3
#


#######################################################################

.is_mem_matrix <- function(m)
{
    is.matrix(m) || is(m, "Matrix")
}

.x_as_Matrix <- function(from, to)
{
    k <- x_type(from)
    if (k == 2L)
    {
        rv <- blockApply(from, function(x) as(x, to),
            grid=rowAutoGrid(from), as.sparse=NA, BPPARAM=NULL)
        do.call(rbind, rv)
    } else {
        rv <- blockApply(from, function(x) as(x, to),
            grid=colAutoGrid(from), as.sparse=NA, BPPARAM=NULL)
        do.call(cbind, rv)
    }
}

x_as_Matrix <- function(from, to)
{
    x_check(from, "Calling SCArray:::x_as_Matrix() with %s ...")
    if (to %in% "Matrix")
    {
        to <- ifelse(is_sparse(from), "sparseMatrix", "matrix")
    }
    .x_as_Matrix(from, to)
}

x_as_sparseMatrix <- function(from, to)
{
    x_check(from, sprintf(
        "Calling SCArray:::x_as_sparseMatrix() with %s (to %s) ...", "%s", to))
    .x_as_Matrix(from, to)
}

x_as_denseMatrix <- function(from, to)
{
    x_check(from, sprintf(
        "Calling SCArray:::x_as_denseMatrix() with %s (to %s) ...", "%s", to))
    k <- x_type(from)
    if (k == 2L)
    {
        rv <- blockApply(from, identity,
            grid=rowAutoGrid(from), as.sparse=NA, BPPARAM=NULL)
        rv <- do.call(rbind, rv)
    } else {
        rv <- blockApply(from, identity,
            grid=colAutoGrid(from), as.sparse=NA, BPPARAM=NULL)
        rv <- do.call(cbind, rv)
    }
    as(rv, to)
}

# setAs(SMatrix, "sparseMatrix", x_as_sparseMatrix)
# setAs(SMatrix, "CsparseMatrix", x_as_sparseMatrix)
# setAs(SMatrix, "RsparseMatrix", x_as_sparseMatrix)
# setAs(SMatrix, "TsparseMatrix", x_as_sparseMatrix)
# setAs(SMatrix, "dgCMatrix", x_as_sparseMatrix)
# setAs(SMatrix, "dgRMatrix", x_as_sparseMatrix)
# setAs(SMatrix, "dgTMatrix", x_as_sparseMatrix)
# setAs(SMatrix, "Matrix", x_as_Matrix)
# setAs(SMatrix, "denseMatrix", x_as_denseMatrix)



#######################################################################

.crossprod_x_row_seq <- function(x)
{
    gd <- rowAutoGrid(x)
    pb <- x_progress(gd)
    if (!is.null(pb)) on.exit(close(pb))
    # block processing
    blockReduce(function(bk, v, pb)
    {
        if (is(bk, "SparseArraySeed")) bk <- as(bk, "sparseMatrix")
        if (!is.null(pb))
            setTxtProgressBar(pb, currentBlockId())
        .Call(c_add, v, as.matrix(crossprod(bk)))
    }, x, matrix(0, ncol(x), ncol(x)), grid=gd, as.sparse=NA, pb=pb)
}


.crossprod_x_row_parallel <- function(x, sp, bp)
{
    # run
    lst <- bplapply(sp, function(s, x)
    {
        dm <- dim(x)
        if (s[1L]!=1L || s[2L]!=dm[1L])
        {
            # sub rows
            x <- x[seq.int(s[1L], s[2L]), , drop=FALSE]
        }
        gd <- rowAutoGrid(x)
        if (s[3L]==1L && s[4L]==dm[2L])
        {
            # crossprod the whole matrix
            blockReduce(function(bk, v, ii)
            {
                if (is(bk, "SparseArraySeed"))
                    bk <- as(bk, "sparseMatrix")
                .Call(c_add, v, as.matrix(crossprod(bk)))
            }, x, matrix(0, ncol(x), ncol(x)), grid=gd, as.sparse=NA)
        } else {
            # crossprod the sub matrix
            ii <- seq.int(s[3L], s[4L])  # sub columns
            blockReduce(function(bk, v, ii)
            {
                if (is(bk, "SparseArraySeed"))
                    bk <- as(bk, "sparseMatrix")
                .Call(c_add, v, as.matrix(crossprod(bk, bk[,ii])))
            }, x, matrix(0, ncol(x), length(ii)), grid=gd, as.sparse=NA, ii=ii)
        }
    }, BPPARAM=bp, x=x)
    # return(list(lst=lst, sp=sp))
    # reduce
    irow <- vapply(sp, `[`, i=3L, 0L)
    iflag <- rep(TRUE, length(irow))
    for (i in unique(irow))
    {
        ii <- which(irow == i)
        k <- ii[1L]
        for (j in ii[-1L])
        {
            lst[[k]] <- .Call(c_add, lst[[k]], lst[[j]])
            iflag[j] <- FALSE
        }
    }
    lst <- lst[iflag]
    # output
    base::Reduce(cbind, lst)
}


# t(x) %*% y, where y = x
x_crossprod_x <- function(x, y)
{
    x_check(x, "Calling SCArray:::x_crossprod_x() with %s ...")
    k <- x_type(x)
    if (k == 2L)
    {
        # efficient direction of 'x' is row, result in [ ncol(x), ncol(x) ]
        bp <- getAutoBPPARAM()
        nworker <- .get_num_worker(bp)
        # split rows & columns (not larger than 2GB)
        ngrp <- floor(min(getAutoBlockSize(), 2^31) / (ncol(x)^2*8))
        sp_r <- scNumSplit(nrow(x), min(ngrp, nworker))
        sp_c <- scNumSplit(ncol(x), floor(nworker/length(sp_r)))
        # do
        if (length(sp_r)!=1L || length(sp_c)!=1L)
        {
            x_msg(sprintf(
                "Distributed to %d processes with (%d x %d) sub matrices ...",
                nworker, length(sp_r), length(sp_c)))
            # split list
            sp <- list()
            for (i in sp_r)
                for (j in sp_c) sp <- c(sp, list(c(i, j)))
            # distribute
            rv <- .crossprod_x_row_parallel(x, sp, bp)
        } else {
            # sequentially
            rv <- .crossprod_x_row_seq(x)
        }
        # output
        s <- colnames(x)
        if (!is.null(s)) dimnames(rv) <- list(s, s)
        rv
    } else {
        crossprod(as(x, DMatrix))
    }
}


################

# x %*% t(y), where y = x
x_tcrossprod_x <- function(x, y)
{
    x_check(x, "Calling SCArray:::x_tcrossprod_x() with %s ...")
    k <- x_type(x)
    if (k == 1L)
    {
        # efficient direction of 'x' is column, result in [ nrow(x), nrow(x) ]
        rv <- blockReduce(function(bk, v)
        {
            if (is(bk, "SparseArraySeed")) bk <- as(bk, "sparseMatrix")
            .Call(c_add, v, as.matrix(tcrossprod(bk)))
        }, x, matrix(0.0, nrow(x), nrow(x)), grid=colAutoGrid(x), as.sparse=NA)
        # output
        s <- rownames(x)
        if (!is.null(s)) dimnames(rv) <- list(s, s)
        rv
    } else {
        tcrossprod(as(x, DMatrix))
    }
}


setMethod("crossprod", c(SMatrix, "missing"), x_crossprod_x)
setMethod("crossprod", c(SMatrix, "ANY"), function(x, y) t(x) %*% y)
setMethod("crossprod", c("ANY", SMatrix), function(x, y) t(x) %*% y)

setMethod("tcrossprod", c(SMatrix, "missing"), x_tcrossprod_x)
setMethod("tcrossprod", c(SMatrix, "ANY"), function(x, y) x %*% t(y))
setMethod("tcrossprod", c("ANY", SMatrix), function(x, y) x %*% t(y))


################

# x %*% y, where x is GDS and y is matrix
x_mul_x_y0 <- function(x, y)
{
    x_msg("\\=> Calling SCArray:::x_mul_x_y0()")
    # initialize
    rnames <- rownames(x)
    cnames <- colnames(y)
    if (!is.null(dimnames(x))) dimnames(x) <- NULL
    # if (!is.null(dimnames(y))) dimnames(y) <- list(NULL, NULL)
    # block processing
    gd <- colAutoGrid(x)
    pb <- x_progress(gd)
    if (!is.null(pb)) on.exit(close(pb))
    rv <- blockReduce(function(bk, v, ym, pb)
    {
        if (is(bk, "SparseArraySeed")) bk <- as(bk, "sparseMatrix")
        vw <- currentViewport()
        ii <- start(vw)[2L]:end(vw)[2L]
        if (!is.null(pb))
            setTxtProgressBar(pb, currentBlockId())
        .Call(c_add, v, as.matrix(bk %*% ym[ii, , drop=FALSE]))
    }, x, matrix(0, nrow(x), ncol(y)), grid=gd, as.sparse=NA, ym=y, pb=pb)
    # output
    rownames(rv) <- rnames
    colnames(rv) <- cnames
    rv
}


x_multiply_x_yANY <- function(x, y)
{
    # check
    x_check(x, "Calling SCArray:::x_multiply_x_yANY() with %s ...")
    if (!.is_mem_matrix(y))
    {
        if (!is.vector(y))
        {
            stop(paste("Matrix multiplication of a", class(x), "by a",
                class(y), "object is not supported"))
        }
        y <- matrix(y, ncol=1L)
    }
    if (ncol(x) != nrow(y))
        stop("non-conformable arguments")
    # compute ...
    k <- x_type(x)
    if (k == 1L)
    {
        # efficient direction of 'x' is column
        x_mul_x_y0(x, y)
    } else {
        as.matrix(as(x, DMatrix) %*% y)
    }
}

x_multiply_xANY_y <- function(x, y)
{
    # check
    x_check(y, "Calling SCArray:::x_multiply_xANY_y() with %s ...")
    if (!.is_mem_matrix(x))
    {
        if (!is.vector(x))
        {
            stop(paste("Matrix multiplication of a", class(x), "object by a",
                class(y), "is not supported"))
        }
        x <- matrix(x, nrow=1L)
    }
    if (ncol(x) != nrow(y))
        stop("non-conformable arguments")
    # compute ...
    k <- x_type(y)
    if (k == 2L)
    {
        # efficient direction of 'y' is row
        t(x_mul_x_y0(t(y), t(x)))
    } else {
        as.matrix(x %*% as(y, DMatrix))
    }
}

x_multiply_x_y <- function(x, y)
{
    as.matrix(as(x, DMatrix) %*% as(y, DMatrix))
}

setMethod("%*%", c(SMatrix, "ANY"), x_multiply_x_yANY)
setMethod("%*%", c("ANY", SMatrix), x_multiply_xANY_y)
# setMethod("%*%", c(SMatrix, SMatrix), x_multiply_x_y)

