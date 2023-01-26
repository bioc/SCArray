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

# t(x) %*% y, y = x
x_crossprod_x <- function(x, y)
{
    x_check(x, "Calling SCArray:::x_crossprod_x() with %s ...")
    k <- x_type(x)
    if (k == 2L)
    {
        # efficient direction of 'x' is Y, result in [ ncol(x), ncol(x) ]
        blockReduce(function(bk, v) {
            if (is(bk, "SparseArraySeed")) bk <- as(bk, "CsparseMatrix")
            crossprod(bk) + v
        }, x, matrix(0.0, nrow=ncol(x), ncol=ncol(x)),
            grid=rowAutoGrid(x), as.sparse=NA)
    } else
        crossprod(as(x, DMatrix))
}

# x %*% t(y), y = x
x_tcrossprod_x <- function(x, y)
{
    x_check(x, "Calling SCArray:::x_tcrossprod_x() with %s ...")
    k <- x_type(x)
    if (k == 1L)
    {
        # efficient direction of 'x' is X, result in [ nrow(x), nrow(x) ]
        blockReduce(function(bk, v) {
            if (is(bk, "SparseArraySeed")) bk <- as(bk, "CsparseMatrix")
            .Call(c_add_update, v, tcrossprod(bk))
        }, x, matrix(0.0, nrow=nrow(x), ncol=nrow(x)),
            grid=colAutoGrid(x), as.sparse=NA)
    } else
        tcrossprod(as(x, DMatrix))
}


setMethod("crossprod", c("SC_GDSMatrix", "missing"), x_crossprod_x)
setMethod("tcrossprod", c("SC_GDSMatrix", "missing"), x_tcrossprod_x)

