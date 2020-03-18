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

# GDSArray and GDSMatrix classes
setClass("GDS_scArray", contains="DelayedArray")
setClass("GDS_scMatrix", contains=c("DelayedMatrix", "GDS_scArray"))

setAs("GDS_scArray", "GDS_scMatrix", function(from) new("GDS_scMatrix", from))
setAs("GDS_scMatrix", "GDS_scArray", function(from) from)
setAs("ANY", "GDS_scMatrix",
    function(from) as(as(from, "GDS_scArray"), "GDS_scMatrix"))

setValidity("GDS_scArray", function(object)
    {
        if (!is(object@seed, "SCArraySeed"))
            return("'object@seed' must be a SCArraySeed object")
        TRUE
    }
)

setMethod("DelayedArray", "SCArraySeed",
    function(seed) new_DelayedArray(seed, Class="GDS_scArray")
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
.extract_array_from_SCArraySeed <- function(x, index)
{
    ans_dim <- DelayedArray:::get_Nindex_lengths(index, dim(x))
    if (any(ans_dim == 0L))
    {
        tp <- as.character(objdesp.gdsn(index.gdsn(x@gds, x@varname))$type)
        ans <- switch(tp,
            Raw=raw(), Integer=integer(), Logical=logical(),
            Real=double(), String=character(),
            stop("Unsupported data type: ", tp))
        dim(ans) <- ans_dim
    } else {
        ans <- readex.gdsn(index.gdsn(x@gds, x@varname), index, .sparse=FALSE)
        if (!is.array(ans))  ## 'ans' must be an array
            dim(ans) <- ans_dim
    }
    ans
}

# GDSArray constructor and coercion methods.
setMethod("extract_array", "SCArraySeed", .extract_array_from_SCArraySeed)


#######################################################################
# SCArraySeed constructor

SCArraySeed <- function(gds, varname)
{
    # check gds
    if (is.character(gds))
    {
        stopifnot(length(gds) == 1L)
        gds <- scOpen(gds, allow.duplicate=TRUE)
        on.exit(scClose(gds))
    }
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
    on.exit()

    # output
    new2("SCArraySeed", gds=gds, varname=varname, dim=dm,
        dimnames=vector("list", length(dm)))
}


scArray <- function(gds, varname)
{
    seed <- SCArraySeed(gds, varname)
    DelayedArray(seed)
}
