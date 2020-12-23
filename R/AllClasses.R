#######################################################################
#
# Package name: SCArray
#
# Description:
#     Large-scale single-cell RNA-seq data manipulation with GDS files
#
# Copyright (C) 2020    Xiuwen Zheng / AbbVie-ComputationalGenomics
# License: GPL-3
#


#######################################################################
# Class definition

setClass("SCArrayFileClass", contains="gds.class")

# test the validity of objects
setValidity("SCArrayFileClass", function(object)
    {
        if (!inherits(object, "gds.class"))
            return("object should inherited from 'gds.class'.")
        var.names <- ls.gdsn(object)
        if (!all(c("feature.id", "sample.id") %in% var.names))
            return("feature.id and sample.id are required variables.")
        TRUE
    }
)


#######################################################################
# Class definition
# See the vignette of the DelayedArray package
#

# if the file is open, no action internally
.reopen <- function(x) gdsfmt:::.reopen(x@gds)

# slot member access
.gds      <- function(x) x@gds
.filename <- function(x) x@filename
.varname  <- function(x) x@varname
.dim      <- function(x) x@dim


# define SCArraySeed class
setClass("SCArraySeed", contains="Array",
    slots = c(
        gds = "SCArrayFileClass",
        filename = "character",
        varname = "character",
        dim  = "integer",
        dimnames = "list"
    )
)


# set the DelayedArray function
setMethod("DelayedArray", "SCArraySeed",
    function(seed) new_DelayedArray(seed, Class="DelayedArray")
)


# show method for SCArraySeed object
setMethod("show", "SCArraySeed", function(object)
    {
        .cat("SCArraySeed\n",
            "File: ", .filename(object), "\n",
            "Array node: ", .varname(object), "\n",
            "Dim: ", paste(.dim(object), collapse=" x "))
    }
)


# extract an array from DelayedArray
setMethod("extract_array", "SCArraySeed", function(x, index)
    {
        # check
        stopifnot(is.list(index), !anyNA(index))
        ans_dim <- DelayedArray:::get_Nindex_lengths(index, dim(x))
        # reopen the file if needed
        .reopen(x)
        # read
        if (any(ans_dim == 0L))
        {
            tp <- objdesp.gdsn(index.gdsn(.gds(x), .varname(x)))$type
            ans <- switch(as.character(tp),
                Raw=raw(), Integer=integer(), Logical=logical(),
                Real=double(), String=character(),
                stop("Unsupported data type: ", tp))
            dim(ans) <- ans_dim
        } else {
            nd <- index.gdsn(.gds(x), .varname(x))
            ans <- readex.gdsn(nd, index, .sparse=FALSE)
            if (!is.array(ans))  # ans must be an array
                dim(ans) <- ans_dim
        }
        ans
    }
)


# return whether the array is sparse or not
setMethod("is_sparse", "SCArraySeed", function(x)
    {
        # reopen the file if needed
        .reopen(x)
        # get data type
        tr <- objdesp.gdsn(index.gdsn(x@gds, x@varname))$trait
        tr %in% c("SparseReal32", "SparseReal64", "SparseInt8", "SparseUInt8",
            "SparseInt16", "SparseUInt16", "SparseInt32", "SparseUInt32",
            "SparseInt64", "SparseUInt64")
    }
)


# extract an array from DelayedArray
setMethod("extract_sparse_array", "SCArraySeed", function(x, index)
    {
        # check
        stopifnot(is.list(index), !anyNA(index))
        ans_dim <- DelayedArray:::get_Nindex_lengths(index, dim(x))
        # reopen the file if needed
        .reopen(x)
        # read
        if (any(ans_dim == 0L))
        {
            tp <- objdesp.gdsn(index.gdsn(.gds(x), .varname(x)))$type
            ans <- switch(as.character(tp),
                Raw=raw(), Integer=integer(), Logical=logical(),
                Real=double(), String=character(),
                stop("Unsupported data type: ", tp))
            SparseArraySeed(ans_dim, nzdata=ans)
        } else {
            nd <- index.gdsn(.gds(x), .varname(x))
            ans <- readex.gdsn(nd, index, .sparse=TRUE)
            if (inherits(ans, "gds_sparse_nz_class"))
            {
                # ans is a list(nzdata, nzindex)
                SparseArraySeed(ans_dim, ans$nzindex, ans$nzdata, check=FALSE)
            } else {
                # ans could be a dense array, dgCMatrix or lgCMatrix
                as(ans, "SparseArraySeed")
            }
        }
    }
)


# SCArraySeed constructor
SCArraySeed <- function(gds, varname)
{
    # check gds
    stopifnot(inherits(gds, "SCArrayFileClass"))
    # check varname
    stopifnot(is.character(varname), length(varname)==1L, !is.na(varname))
    nd <- index.gdsn(gds, varname, silent=TRUE)
    if (is.null(nd))
        stop("No '", varname, "'")
    # check dimension
    dp <- objdesp.gdsn(nd)
    if (!dp$is.array)
        stop("'", varname, "' is not an array.")
    dm <- dp$dim
    # output
    new2("SCArraySeed", gds=gds, filename=gds$filename, varname=varname,
        dim=dm, dimnames=vector("list", length(dm)))
}
