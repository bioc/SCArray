#######################################################################
#
# Package name: SCArray
#
# Description:
#     Large-scale single-cell RNA-seq data manipulation with GDS files
#
# Copyright (C) 2021    Xiuwen Zheng (@AbbVie-ComputationalGenomics)
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
setAs(
    "ANY", "SC_GDSMatrix",
    function(from) as(as(from, "SC_GDSArray"), "SC_GDSMatrix"))

# setMethod(
#     "DelayedArray", "SCArraySeed",
#     function(seed) new_DelayedArray(seed, Class="SC_GDSArray")
# )



#######################################################################

.sc_rowVars <- function(x, rows=NULL, cols=NULL, na.rm=FALSE, ...)
{
    # check
    DelayedArray:::.get_ans_type(x, must.be.numeric=TRUE)

    # block read
    pm <- blockReduce(function(bk, v) {
        .Call(c_rowVars_update, bk, v)
    }, x, init=rep(0.0, nrow(x)*3L), grid=colAutoGrid(x))

    # output
    .Call(c_rowVars_final, pm)
}

setMethod("rowVars", "SC_GDSMatrix", .sc_rowVars)













