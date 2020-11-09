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
scOpen <- function(gdsfn, readonly=TRUE, allow.duplicate=TRUE)
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
    if (is.character(gdsfile))
        gdsfile <- scOpen(gdsfile, readonly=TRUE, allow.duplicate=TRUE)
    stopifnot(inherits(gdsfile, "SCArrayFileClass"))
    # new DelayedArray
    seed <- SCArraySeed(gdsfile, varname)
    DelayedArray(seed)
}


#######################################################################
# Get an SingleCellExperiment/SummarizedExperiment instance
#
scExperiment <- function(gdsfile, sce=TRUE, use.names=TRUE, load.row=TRUE,
    load.col=TRUE)
{
    # check
    if (is.character(gdsfile))
        gdsfile <- scOpen(gdsfile, readonly=TRUE, allow.duplicate=TRUE)
    stopifnot(inherits(gdsfile, "SCArrayFileClass"))
    stopifnot(is.logical(sce), length(sce)==1L)
    stopifnot(is.logical(use.names), length(use.names)==1L)
    stopifnot(is.logical(load.row), length(load.row)==1L)
    stopifnot(is.logical(load.col), length(load.col)==1L)

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

    # rowData
    rowdat <- NULL
    if (isTRUE(load.row))
    {
        nd <- index.gdsn(gdsfile, "feature.data", silent=TRUE)
        if (!is.null(nd))
        {
            nmlst <- ls.gdsn(nd)
            v <- lapply(nmlst, function(nm) read.gdsn(index.gdsn(nd, nm)))
            names(v) <- nmlst
            rowdat <- DataFrame(v, row.names=feat_id)
        }
    }

    # colData
    coldat <- NULL
    if (isTRUE(load.col))
    {
        nd <- index.gdsn(gdsfile, "sample.data", silent=TRUE)
        if (!is.null(nd))
        {
            nmlst <- ls.gdsn(nd)
            v <- lapply(nmlst, function(nm) read.gdsn(index.gdsn(nd, nm)))
            names(v) <- nmlst
            coldat <- DataFrame(v, row.names=samp_id)
        }
    }

    # output
    if (isTRUE(sce))
    {
        if (is.null(coldat))
            SingleCellExperiment(assays=lst, rowData=rowdat)
        else
            SingleCellExperiment(assays=lst, rowData=rowdat, colData=coldat)
    } else {
        if (is.null(coldat))
            SummarizedExperiment(assays=lst, rowData=rowdat)
        else
            SummarizedExperiment(assays=lst, rowData=rowdat, colData=coldat)
    }
}


#######################################################################
# Convert an R object (count matrix) to a single-cell GDS file
# the R object can be matrix, DelayedMatrix
#
scConvGDS <- function(obj, outfn, save.sp=TRUE,
    type=c("float32", "float64", "int32"), compress="LZMA_RA", clean=TRUE,
    verbose=TRUE)
{
    # check
    stopifnot(is.matrix(obj) | is(obj, "DelayedMatrix") |
        is(obj, "SingleCellExperiment") | is(obj, "SummarizedExperiment"))
    stopifnot(is.character(outfn), length(outfn)==1L)
    stopifnot(is.logical(save.sp), length(save.sp)==1L)
    stopifnot(is.character(compress), length(compress)==1L)
    stopifnot(is.logical(clean), length(clean)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)
    type <- match.arg(type)

    # row and column names
    rownm <- colnm <- rowdat <- coldat <- metadat <- NULL
    nr <- nc <- 0L
    lst <- list()
    if (is.matrix(obj) || is(obj, "DelayedMatrix"))
    {
        rownm <- rownames(obj); colnm <- colnames(obj)
        nr <- nrow(obj); nc <- ncol(obj)
        assaylst$counts <- obj
    } else if (is(obj, "SingleCellExperiment") |
        is(obj, "SummarizedExperiment"))
    {
        rownm <- rownames(obj)
        if (is.null(rownm)) stop("No rownames!")
        colnm <- colnames(obj)
        if (is.null(colnm)) stop("No colnames!")
        nr <- length(rownm); nc <- length(colnm)
        assaylst <- assays(obj)
        rowdat <- rowData(obj)
        coldat <- colData(obj)
        metadat <- metadata(obj)
    }

    if (is.null(rownm)) rownm <- seq_len(nr)
    if (is.null(colnm)) colnm <- seq_len(nc)

    # create a gds file
    if (verbose) .cat("Output: ", outfn)
    outf <- createfn.gds(outfn)
    on.exit(closefn.gds(outf))
    put.attr.gdsn(outf$root, "FileFormat", "SC_ARRAY")
    put.attr.gdsn(outf$root, "FileVersion", "v1.0")
    if (verbose) .cat("Compression: ", compress)

    add.gdsn(outf, "feature.id", rownm, compress=compress, closezip=TRUE)
    add.gdsn(outf, "sample.id", colnm, compress=compress, closezip=TRUE)
    if (verbose) .cat("Dimension: ", nr, " x ", nc)

    if (verbose) .cat("Assay List:")
    for (i in seq_along(assaylst))
    {
        nm <- names(assaylst)[i]
        if (verbose) cat("   ", nm)
        st <- type
        if (save.sp)
        {
            st <- switch(type, float32="sp.real32", float64="sp.real64",
                int32="sp.int32", stop("Invalid 'type': ", type))
        }
        nd <- add.gdsn(outf, nm, valdim=c(nr, 0L), compress=compress,
            storage=st)
        mt <- assaylst[[i]]
        blockApply(mt, function(x) { append.gdsn(nd, x); NULL },
            grid=colAutoGrid(mt), BPPARAM=NULL)
        readmode.gdsn(nd)
        if (verbose) { cat("  |"); print(nd) }
    }

    nfd <- addfolder.gdsn(outf, "feature.data")
    if (!is.null(rowdat))
    {
        if (verbose) .cat("rowData:")
        for (i in seq_len(ncol(rowdat)))
        {
            nm <- names(rowdat)[i]
            .cat("    ", nm)
            add.gdsn(nfd, nm, rowdat[,i], compress=compress, closezip=TRUE)
        }
    }

    nfd <- addfolder.gdsn(outf, "sample.data")
    if (!is.null(coldat))
    {
        if (verbose) .cat("colData:")
        for (i in seq_len(ncol(coldat)))
        {
            nm <- names(coldat)[i]
            nm <- gsub("/", "_", nm, fixed=TRUE)
            .cat("    ", nm)
            add.gdsn(nfd, nm, coldat[,i], compress=compress, closezip=TRUE)
        }
    }

    nfd <- addfolder.gdsn(outf, "meta.data")
    if (length(metadat) > 0L)
    {
        if (verbose) .cat("metadata:")
        for (i in seq_along(metadat))
        {
            nm <- names(metadat)[i]
            .cat("    ", nm)
            v <- metadat[[i]]
            if (is(v, "DataFrame")) v <- as.data.frame(v)
            add.gdsn(nfd, nm, v, compress=compress, closezip=TRUE)
        }
    }

    # close file
    on.exit()
    closefn.gds(outf)
    .cat("Done.")
    if (isTRUE(clean)) cleanup.gds(outfn)

    # output
    invisible(normalizePath(outfn))
}
