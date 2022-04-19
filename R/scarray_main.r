#######################################################################
#
# Package name: SCArray
#
# Description:
#     Large-scale single-cell RNA-seq data manipulation with GDS files
#
# Copyright (C) 2020-2022    Xiuwen Zheng (@AbbVie-ComputationalGenomics)
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

    # check the file format
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

    # output
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
    x <- vapply(nm, FUN=function(s) {
        identical(objdesp.gdsn(index.gdsn(gdsfile, s))$dim, dm)
    }, FUN.VALUE=TRUE)
    lst <- lapply(nm[x], function(s) {
        m <- scArray(gdsfile, s)
        rownames(m) <- feat_id; colnames(m) <- samp_id
        if (length(dim(m))==2L)
            m <- as(m, "SC_GDSMatrix")
        else
            m <- as(m, "SC_GDSArray")
        m
    })
    names(lst) <- nm[x]

    # load rowData
    rowdat <- NULL
    if (isTRUE(load.row))
    {
        nd <- index.gdsn(gdsfile, "feature.data", silent=TRUE)
        if (!is.null(nd))
        {
            nmlst <- ls.gdsn(nd, include.hidden=FALSE)
            v <- lapply(nmlst, function(nm) read.gdsn(index.gdsn(nd, nm)))
            names(v) <- nmlst
            rowdat <- DataFrame(v, row.names=feat_id)
        }
    }

    # load colData
    coldat <- NULL
    if (isTRUE(load.col))
    {
        nd <- index.gdsn(gdsfile, "sample.data", silent=TRUE)
        if (!is.null(nd))
        {
            nmlst <- ls.gdsn(nd, include.hidden=FALSE)
            v <- lapply(nmlst, function(nm) read.gdsn(index.gdsn(nd, nm)))
            names(v) <- nmlst
            coldat <- DataFrame(v, row.names=samp_id)
        }
    }

    # output
    if (isTRUE(sce))
    {
        # return a SingleCellExperiment object
        if (is.null(coldat))
            SingleCellExperiment(assays=lst, rowData=rowdat)
        else
            SingleCellExperiment(assays=lst, rowData=rowdat, colData=coldat)
    } else {
        # return a SummarizedExperiment object
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
    stopifnot(is.matrix(obj) | inherits(obj, "Matrix") |
        is(obj, "DelayedMatrix") | is(obj, "SingleCellExperiment") |
        is(obj, "SummarizedExperiment"))
    stopifnot(is.character(outfn), length(outfn)==1L)
    stopifnot(is.logical(save.sp), length(save.sp)==1L)
    type <- match.arg(type)
    stopifnot(is.character(compress), length(compress)==1L)
    stopifnot(is.logical(clean), length(clean)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    # row and column names
    rownm <- colnm <- rowdat <- coldat <- metadat <- NULL
    nr <- nc <- 0L
    lst <- list()
    if (is.matrix(obj) || inherits(obj, "Matrix") || is(obj, "DelayedMatrix"))
    {
        rownm <- rownames(obj); colnm <- colnames(obj)
        nr <- nrow(obj); nc <- ncol(obj)
        assaylst <- list(counts=obj)
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

    # need rownames and colnames
    if (is.null(rownm)) rownm <- seq_len(nr)
    if (anyDuplicated(rownm)) stop("rownames should be unique.")
    if (is.null(colnm)) colnm <- seq_len(nc)
    if (anyDuplicated(colnm)) stop("colnames should be unique.")

    # create a gds file
    if (verbose) .cat("Output: ", outfn)
    outf <- createfn.gds(outfn)
    on.exit(closefn.gds(outf))
    put.attr.gdsn(outf$root, "FileFormat", "SC_ARRAY")
    put.attr.gdsn(outf$root, "FileVersion", "v1.0")
    if (verbose) .cat("Compression: ", compress)

    # add feature and sample IDs
    add.gdsn(outf, "feature.id", rownm, compress=compress, closezip=TRUE)
    add.gdsn(outf, "sample.id", colnm, compress=compress, closezip=TRUE)
    if (verbose) .cat("Dimension: ", nr, " x ", nc)

    # assays
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

    # add rowData to feature.data
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

    # add colData to sample.data
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

    # add metadata to meta.data
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
    if (isTRUE(clean)) cleanup.gds(outfn, verbose=FALSE)
    cat("Done.\n")

    # output
    invisible(normalizePath(outfn))
}


#######################################################################
# Convert CellRanger Market Exchange Format (MEX) files to GDS
#
scMEX2GDS <- function(feature_fn, barcode_fn, mtx_fn, outfn,
    feature_colnm=c("id", "gene", "feature_type"),
    type=c("float32", "float64", "int32"), compress="LZMA_RA", clean=TRUE,
    verbose=TRUE)
{
    # check
    stopifnot(is.character(feature_fn), length(feature_fn)==1L)
    stopifnot(is.character(barcode_fn), length(barcode_fn)==1L)
    stopifnot(is.character(mtx_fn), length(mtx_fn)==1L)
    stopifnot(is.character(outfn), length(outfn)==1L)
    stopifnot(is.character(feature_colnm))
    type <- match.arg(type)
    stopifnot(is.character(compress), length(compress)==1L)
    stopifnot(is.logical(clean), length(clean)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    # load features
    if (anyDuplicated(feature_colnm))
        stop("'feature_colnm' should be unique.")
    if (length(feature_colnm) == 0L)
        feature_colnm <- "id"
    if (verbose) .cat("Load ", sQuote(feature_fn))
    ft <- read.delim(feature_fn, header=FALSE, stringsAsFactors=FALSE)
    if (ncol(ft) > length(feature_colnm))
    {
        feature_colnm <- c(feature_colnm,
            paste0("var", seq.int(length(feature_colnm)+1L,
                length.out=ncol(ft)-length(feature_colnm))))
    }
    colnames(ft) <- feature_colnm[seq_len(ncol(ft))]
    if (!("id" %in% colnames(ft)))
        stop("No id found in the feature file.")
    if (anyDuplicated(ft$id))
        stop("The first column of feature file is used as feature.id, and should be unique.")

    # load barcodes
    if (verbose) .cat("Load ", sQuote(barcode_fn))
    bt <- read.delim(barcode_fn, header=FALSE, stringsAsFactors=FALSE)
    if (anyDuplicated(bt$V1))
        stop("The first column of barcode file is used as sample.id, and should be unique.")

    # load count matrix
    if (verbose) .cat("Load ", sQuote(mtx_fn))
    mt <- Matrix::readMM(mtx_fn)
    if (verbose)
    {
        nz <- Matrix::nnzero(mt)
        .cat("    ", nrow(mt), "x", ncol(mt), ", # of nonzeros: ",
            nz, sprintf(" (%.4f%%)", 100*nz/prod(dim(mt))))
    }
    if (nrow(mt) != nrow(ft))
        stop("# of rows should be as the same as feature file.")
    if (ncol(mt) != nrow(bt))
        stop("# of columns should be as the same as # of rows in barcode file.")

    # create a gds file
    if (verbose) .cat("Output: ", outfn)
    outf <- createfn.gds(outfn)
    on.exit(closefn.gds(outf))
    put.attr.gdsn(outf$root, "FileFormat", "SC_ARRAY")
    put.attr.gdsn(outf$root, "FileVersion", "v1.0")
    if (verbose) .cat("Compression: ", compress)

    # add feature and sample IDs
    add.gdsn(outf, "feature.id", ft$id, compress=compress, closezip=TRUE)
    add.gdsn(outf, "sample.id", bt$V1, compress=compress, closezip=TRUE)

    # assays
    if (verbose) .cat("Count matrix ...")
    st <- switch(type, float32="sp.real32", float64="sp.real64",
        int32="sp.int32", stop("Invalid 'type': ", type))
    nd <- add.gdsn(outf, "counts", valdim=c(nrow(mt), 0L), compress=compress,
        storage=st)
    blockApply(mt, function(x) { append.gdsn(nd, x); NULL },
        grid=colAutoGrid(mt), BPPARAM=NULL)
    readmode.gdsn(nd)
    if (verbose) { cat("  |"); print(nd) }

    # add rowData to feature.data
    nfd <- addfolder.gdsn(outf, "feature.data")
    for (i in seq.int(2L, length.out=ncol(ft)-1L))
    {
        add.gdsn(nfd, names(ft)[i], ft[[i]], compress=compress,
            closezip=TRUE)
    }

    # add colData to sample.data
    nfd <- addfolder.gdsn(outf, "sample.data")
    for (i in seq_len(ncol(bt)-1L))
    {
        add.gdsn(nfd, paste0("var", i), bt[[i+1L]], compress=compress,
            closezip=TRUE)
    }

    # close file
    on.exit()
    closefn.gds(outf)
    if (isTRUE(clean)) cleanup.gds(outfn, verbose=FALSE)
    cat("Done.\n")

    # output
    invisible(normalizePath(outfn))
}


#######################################################################
# Convert CellRanger HDF5 files to GDS
#
scHDF2GDS <- function(h5_fn, outfn, group="matrix", feature_path=character(),
    type=c("float32", "float64", "int32"), compress="LZMA_RA", clean=TRUE,
    verbose=TRUE)
{
    # check
    stopifnot(is.character(h5_fn), length(h5_fn)==1L)
    stopifnot(is.character(outfn), length(outfn)==1L)
    stopifnot(is.character(group), length(group)==1L)
    stopifnot(is.character(feature_path))
    if (anyDuplicated(feature_path))
        stop("'feature_path' should be unique.")
    type <- match.arg(type)
    stopifnot(is.character(compress), length(compress)==1L)
    stopifnot(is.logical(clean), length(clean)==1L)
    stopifnot(is.logical(verbose), length(verbose)==1L)

    # load HDF5 count matrix
    sep <- ifelse(grepl("/$", group), "", "/")
    if (verbose)
    {
        .cat("Load ", sQuote(h5_fn))
        cat("   ", group)
    }
    mt <- HDF5Array::H5SparseMatrix(h5_fn, group)
    if (verbose)
        .cat("  [", class(mt), ": ", nrow(mt), "x", ncol(mt), "]")

    # load barcodes
    nm <- paste(group, "barcodes", sep=sep)
    if (verbose) .cat("    ", nm)
    bt <- rhdf5::h5read(h5_fn, nm)
    if (anyDuplicated(bt))
        stop("'barcodes' should be unique.")
    if (length(bt) != ncol(mt))
        stop("Num. of columns in the count matrix should be # of barcodes.")

    # load features
    if (length(feature_path) == 0L)
    {
        feature_path <- c("genes", "gene_names",
            "features/id", "features/name", "features/feature_type",
            "features/genome")
    }
    a <- rhdf5::h5ls(h5_fn, all=T)
    a$path <- substring(paste(a$group, a$name, sep="/"), 2L)
    test_path <- paste(group, feature_path, sep="/")
    used <- intersect(test_path, a$path)
    if (length(used) == 0L)
        stop("'No valid 'feature_path'.")
    ft <- lapply(used, function(nm) {
        .cat("    ", nm)
        v <- rhdf5::h5read(h5_fn, nm)
        v <- c(v)  # force to regular data type
        if (length(v) != nrow(mt))
            stop(sQuote(nm), " should have ", nrow(mt), " values.")
        v
    })
    names(ft) <- basename(used)
    if (anyDuplicated(ft[[1L]]))
        stop(sQuote(used[1L]), " should be unique.")

    # create a gds file
    if (verbose) .cat("Output: ", outfn)
    outf <- createfn.gds(outfn)
    on.exit(closefn.gds(outf))
    put.attr.gdsn(outf$root, "FileFormat", "SC_ARRAY")
    put.attr.gdsn(outf$root, "FileVersion", "v1.0")
    if (verbose) .cat("Compression: ", compress)

    # add feature and sample IDs
    add.gdsn(outf, "feature.id", ft[[1L]], compress=compress, closezip=TRUE)
    add.gdsn(outf, "sample.id", bt, compress=compress, closezip=TRUE)

    # assays
    if (verbose) .cat("Count matrix ...")
    st <- switch(type, float32="sp.real32", float64="sp.real64",
        int32="sp.int32", stop("Invalid 'type': ", type))
    nd <- add.gdsn(outf, "counts", valdim=c(nrow(mt), 0L), compress=compress,
        storage=st)
    blockApply(mt, function(x) { append.gdsn(nd, x); NULL },
        grid=colAutoGrid(mt), BPPARAM=NULL)
    readmode.gdsn(nd)
    if (verbose) { cat("  |"); print(nd) }

    # add rowData to feature.data
    nfd <- addfolder.gdsn(outf, "feature.data")
    for (i in seq.int(2L, length.out=length(ft)-1L))
    {
        add.gdsn(nfd, names(ft)[i], ft[[i]], compress=compress,
            closezip=TRUE)
    }

    # add colData to sample.data
    addfolder.gdsn(outf, "sample.data")

    # close file
    on.exit()
    closefn.gds(outf)
    if (isTRUE(clean)) cleanup.gds(outfn, verbose=FALSE)
    cat("Done.\n")

    # output
    invisible(normalizePath(outfn))
}


#######################################################################
#
#
#
scObj <- function(obj, verbose=TRUE)
{
    # check
    stopifnot(is.logical(verbose), length(verbose)==1L)

    if (is(obj, "SummarizedExperiment"))
    {
        lst <- assays(obj)
        nm <- names(lst)
        for (i in seq_along(lst))
        {
            v <- lst[[i]]
            if (is(v, "DelayedMatrix") && !is(v, "SC_GDSMatrix"))
            {
                assay(obj, i) <- as(v, "SC_GDSMatrix")
                if (verbose)
                    cat(nm[i], "==> SC_GDSMatrix\n")
            } else if (is(v, "DelayedArray") && !is(v, "SC_GDSArray"))
            {
                assay(obj, i) <- as(v, "SC_GDSArray")
                if (verbose)
                    cat(nm[i], "==> SC_GDSArray\n")
            }
        }
        obj
    } else if (is(obj, "DelayedArray"))
    {
        if (is(obj, "DelayedMatrix") && !is(obj, "SC_GDSMatrix"))
        {
            obj <- as(obj, "SC_GDSMatrix")
            if (verbose)
                cat("==> SC_GDSMatrix\n")
        } else if (is(v, "DelayedArray") && !is(v, "SC_GDSArray"))
        {
            obj <- as(v, "SC_GDSArray")
            if (verbose)
                cat("==> SC_GDSArray\n")
        }
        obj
    } else
        stop("obj should be a SummarizedExperiment object.")
}

