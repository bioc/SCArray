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


# Package-wide variable
.packageEnv <- new.env()


#######################################################################
# Internal functions
#

.cat <- function(...) cat(..., "\n", sep="")

.plural <- function(num) if (num > 1L) "s" else ""

.pretty <- function(x) prettyNum(x, big.mark=",", scientific=FALSE)

.pretty_size <- function(x)
{
    if (x >= 1024^4)
        sprintf("%.1fT", x / 1024^4)
    else if (x >= 1024^3)
        sprintf("%.1fG", x / 1024^3)
    else if (x >= 1024^2)
        sprintf("%.1fM", x / 1024^2)
    else if (x >= 1024L)
        sprintf("%.1fK", x / 1024)
    else if (x > 1L)
        sprintf("%g bytes", x)
    else
        "0 byte"
}


.cfunction <- function(name)
{
    fn <- function(x) NULL
    f <- quote(.Call(SC_ExternalName1, x))
    f[[1L]] <- .Call
    f[[2L]] <- getNativeSymbolInfo(name, "SCArray")$address
    body(fn) <- f
    fn
}

.cfunction2 <- function(name)
{
    fn <- function(x, y) NULL
    f <- quote(.Call(SC_ExternalName2, x, y))
    f[[1L]] <- .Call
    f[[2L]] <- getNativeSymbolInfo(name, "SCArray")$address
    body(fn) <- f
    fn
}



# return total # of features, and total # of samples
.gettotnum <- function(gdsfile)
{
    nfeat <- objdesp.gdsn(index.gdsn(gdsfile, "feature.id"))$dim
    if (length(nfeat) != 1L)
        stop("Invalid dimension of 'feature.id'.")
    nsamp <- objdesp.gdsn(index.gdsn(gdsfile, "sample.id"))$dim
    if (length(nsamp) != 1L)
        stop("Invalid dimension of 'sample.id'.")
    c(nfeat, nsamp)
}



#######################################################################
# Open a SCArray GDS file
#
scOpen <- function(gdsfn, readonly=TRUE, allow.duplicate=FALSE)
{
    # check
    stopifnot(is.character(gdsfn), length(gdsfn)==1L)
    stopifnot(is.logical(readonly), length(readonly)==1L)
    stopifnot(is.logical(allow.duplicate), length(allow.duplicate)==1L)

    # open the file
    ans <- openfn.gds(gdsfn, readonly=readonly, allow.fork=TRUE,
        allow.duplicate=allow.duplicate)

    # FileFormat
    a <- get.attr.gdsn(ans$root)
    if (!is.null(a$FileFormat))
    {
        if (identical(a$FileFormat, "SNP_ARRAY"))
            stop("It is a SNP GDS file, please use SNPRelate::snpgdsOpen().")
        if (identical(a$FileFormat, "SEQ_ARRAY"))
            stop("It is a SeqArray GDS file, please use SeqArray::seqOpen().")
        if (!identical(a$FileFormat, "SC_ARRAY"))
            stop("'FileFormat' should be 'SC_ARRAY'")
    }

    # FileVersion
    version <- a$FileVersion
    if (!is.null(version))
    {
        if (!identical(version, "v1.0"))
        {
            stop(sprintf(
                "FileVersion '%s' (should be v1.0), consider using the updated SCArray",
                as.character(version)))
        }
    }

    # .Call(SEQ_File_Init, ans)
    new("SCArrayFileClass", ans)
}


#######################################################################
# Close the SCArray GDS file
#
scClose <- function(gdsfile)
{
    # check
    stopifnot(inherits(gdsfile, "SCArrayFileClass"))
    # close the GDS file
    closefn.gds(gdsfile)
}


#######################################################################
# Get an DelayedArray instance
#
scArray <- function(gdsfile, varname)
{
    # check
    stopifnot(inherits(gdsfile, "SCArrayFileClass"))
    # new DelayedArray
    seed <- SCArraySeed(gdsfile, varname)
    DelayedArray(seed)
}


#######################################################################
# Get an SingleCellExperiment/SummarizedExperiment instance
#
scExperiment <- function(gdsfile, sc=TRUE, use.names=TRUE)
{
    # check
    stopifnot(inherits(gdsfile, "SCArrayFileClass"))
    stopifnot(is.logical(sc), length(sc)==1L)
    stopifnot(is.logical(use.names), length(use.names)==1L)

    # dimnames
    dm <- .gettotnum(gdsfile)
    feat_id <- samp_id <- NULL
    if (use.names)
    {
        feat_id <- read.gdsn(index.gdsn(gdsfile, "feature.id"))
        samp_id <- read.gdsn(index.gdsn(gdsfile, "sample.id"))
    }
    # list all assays
    nm <- ls.gdsn(gdsfile, include.dirs=FALSE)
    x <- sapply(nm, function(s) {
        dp <- objdesp.gdsn(index.gdsn(gdsfile, s))
        identical(dp$dim, dm)
    })
    lst <- lapply(nm[x], function(s) {
        m <- scArray(gdsfile, s)
        rownames(m) <- feat_id; colnames(m) <- samp_id
        m
    })
    names(lst) <- nm[x]

    # output
    if (isTRUE(sc))
    {
        SingleCellExperiment(assays=lst)
    } else {
        SummarizedExperiment(assays=lst)
    }
}
