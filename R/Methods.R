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

x_runsvd <- function(x, rank, scale=FALSE, approx=TRUE)
{
    # to use crossprod
    if (x_verbose())
        cat("Using cross product for PCA\n")
    # fold = 1 for using cross product
    if (approx)
    {
        rv <- BiocSingular::runIrlbaSVD(x, k=rank, nu=rank, nv=rank,
            center=TRUE, scale=scale, deferred=FALSE, fold=1)
    } else {
        rv <- BiocSingular::runExactSVD(x, k=rank, nu=rank, nv=rank,
            center=TRUE, scale=scale, deferred=FALSE, fold=1)
    }
    # output
    out <- list(sdev = rv$d/sqrt(nrow(x) - 1))
	out$rotation <- rv$v
	colnames(out$rotation) <- sprintf("PC%i", seq_len(ncol(out$rotation)))
    out$x <- sweep(rv$u, 2, rv$d, "*")
    colnames(out$x) <- sprintf("PC%i", seq_len(ncol(out$x)))
	out
}


scRunPCA <- function(x, ncomponents=50, ntop=500, subset_row=NULL, scale=FALSE,
    altexp=NULL, name="PCA", exprs_values="logcounts", dimred=NULL,
    n_dimred=NULL, BSPARAM=NULL, BPPARAM=SerialParam(), verbose=TRUE)
{
    # check
    stopifnot(is(x, "SingleCellExperiment") || is(x, "SummarizedExperiment"))
    stopifnot(is.numeric(ncomponents), length(ncomponents)==1L)

    # get the working matrix 'mat'
    y <- x
    if (!is.null(altexp)) y <- altExp(x, altexp)
    transposed <- !is.null(dimred)
    if (!is.null(dimred))
    {
        mat <- reducedDim(y, dimred)
        if (!is.null(n_dimred))
        {
            if (length(n_dimred) == 1L) n_dimred <- seq_len(n_dimred)
            mat <- mat[, n_dimred, drop = FALSE]
        }
    } else {
        mat <- assay(y, exprs_values)
    }

    # save BPPARAM
    oldbp <- getAutoBPPARAM()
    setAutoBPPARAM(BPPARAM)
    on.exit(setAutoBPPARAM(oldbp))
    if (!bpisup(BPPARAM) && !is(BPPARAM, "MulticoreParam"))
    {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add=TRUE)
    }
    if (!transposed)
    {
        if (verbose)
            .cat("Select top ", ntop, " features with the highest variances ...")
        cv <- rowVars(mat)
        if (is.null(subset_row))
            subset_row <- sort(head(order(cv, decreasing=TRUE), ntop))
        mat <- mat[subset_row, , drop=FALSE]
        cv <- cv[subset_row]
        if (scale)
        {
            flag <- cv >= 1e-08
            mat <- mat[flag, , drop=FALSE]/sqrt(cv[flag])
            cv <- rep(1, nrow(mat))
        }
        mat <- t(mat)
    } else {
        cv <- colVars(mat)
    }

    # do PCA
    if (verbose) cat("Start PCA on the matrix ...\n")
    if (is.null(BSPARAM))
    {
        pca <- x_runsvd(mat, ncomponents, scale=scale)
    } else {
        # use BiocSingular::runPCA instead
        pca <- runPCA(mat, rank=ncomponents, BSPARAM=BSPARAM, BPPARAM=BPPARAM)
    }

    # PCA results
    varExplained <- pca$sdev^2
    percentVar <- varExplained/sum(cv) * 100
    pcs <- pca$x
    rownames(pcs) <- rownames(mat)
    attr(pcs, "varExplained") <- varExplained
    attr(pcs, "percentVar") <- percentVar
    rownames(pca$rotation) <- colnames(mat)
    attr(pcs, "rotation") <- pca$rotation

    # output
    reducedDim(x, name) <- pcs
    x
}

