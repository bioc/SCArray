#######################################################################
#
# Package name: SCArray
#
# Description:
#     Large-scale single-cell RNA-seq data manipulation and analysis with
# efficient algorithms
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

setClass("SCArraySeed", contains="Array",
    slots = c(
        gds = "SCArrayFileClass",
        varname = "character",
        dim  = "integer",
        dimnames = "list"
    )
)

setMethod("DelayedArray", "SCArraySeed",
    function(seed) new_DelayedArray(seed, Class="DelayedArray")
)


# show method for SCArraySeed object
setMethod("show", "SCArraySeed", function(object)
    {
        .cat("SCArraySeed\n",
            "File: ", object@gds$filename, "\n",
            "Array node: ", object@varname, "\n",
            "Dim: ", paste(object@dim, collapse=" x "))
    }
)


# extract an array from DelayedArray
setMethod("extract_array", "SCArraySeed", function(x, index)
    {
        ans_dim <- DelayedArray:::get_Nindex_lengths(index, dim(x))
        if (any(ans_dim == 0L))
        {
            tp <- objdesp.gdsn(index.gdsn(x@gds, x@varname))$type
            ans <- switch(as.character(tp),
                Raw=raw(), Integer=integer(), Logical=logical(),
                Real=double(), String=character(),
                stop("Unsupported data type: ", tp))
            dim(ans) <- ans_dim
        } else {
            nd <- index.gdsn(x@gds, x@varname)
            ans <- readex.gdsn(nd, index, .sparse=FALSE)
            if (!is.array(ans))  ## 'ans' must be an array
                dim(ans) <- ans_dim
        }
        ans
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
    new2("SCArraySeed", gds=gds, varname=varname, dim=dm,
        dimnames=vector("list", length(dm)))
}

