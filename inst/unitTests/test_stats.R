#######################################################################
#
# Consistency regression tests for the SC_GDSMatrix statistical methods.
#
# Each overridden row*/col* method is checked against the matrixStats
# reference computed on the in-memory dense matrix, for both the native
# (type-1) and the transposed (type-2) dispatch paths.
#
#######################################################################

suppressPackageStartupMessages({
    library(RUnit)
    library(SCArray)
    library(DelayedMatrixStats)
    library(SingleCellExperiment)
    library(SummarizedExperiment)
    # NB: do NOT attach matrixStats here -- it would shadow the MatrixGenerics
    # S4 generics with its own matrix-only functions and break dispatch on
    # SC_GDSMatrix. Use matrixStats:: explicitly where a matrix reference is
    # needed.
})

# tolerance (counts are stored as 32-bit float; reductions accumulate in double)
.TOL <- 1e-5


# Compare every supported statistic of an SC_GDSMatrix 'x' against the dense
# reference 'r' (a base matrix holding identical values).
.check_all_stats <- function(x, r, na.rm=FALSE, tag="")
{
    msg <- function(s) paste0(tag, ": ", s)
    # strip dimnames on both sides so these checks compare values only
    # (name handling is verified separately in test_useNames)
    dimnames(x) <- NULL
    dimnames(r) <- NULL

    # -- sums --
    checkEquals(unname(rowSums(x, na.rm=na.rm)),  rowSums(r, na.rm=na.rm),  msg("rowSums"),  tolerance=.TOL)
    checkEquals(unname(colSums(x, na.rm=na.rm)),  colSums(r, na.rm=na.rm),  msg("colSums"),  tolerance=.TOL)
    checkEquals(rowSums2(x, na.rm=na.rm), rowSums2(r, na.rm=na.rm), msg("rowSums2"), tolerance=.TOL)
    checkEquals(colSums2(x, na.rm=na.rm), colSums2(r, na.rm=na.rm), msg("colSums2"), tolerance=.TOL)

    # -- means --
    checkEquals(unname(rowMeans(x, na.rm=na.rm)), rowMeans(r, na.rm=na.rm), msg("rowMeans"), tolerance=.TOL)
    checkEquals(unname(colMeans(x, na.rm=na.rm)), colMeans(r, na.rm=na.rm), msg("colMeans"), tolerance=.TOL)
    checkEquals(rowMeans2(x, na.rm=na.rm), rowMeans2(r, na.rm=na.rm), msg("rowMeans2"), tolerance=.TOL)
    checkEquals(colMeans2(x, na.rm=na.rm), colMeans2(r, na.rm=na.rm), msg("colMeans2"), tolerance=.TOL)

    # -- variance / sd --
    checkEquals(rowVars(x, na.rm=na.rm), rowVars(r, na.rm=na.rm), msg("rowVars"), tolerance=.TOL)
    checkEquals(colVars(x, na.rm=na.rm), colVars(r, na.rm=na.rm), msg("colVars"), tolerance=.TOL)
    checkEquals(rowSds(x, na.rm=na.rm),  rowSds(r, na.rm=na.rm),  msg("rowSds"),  tolerance=.TOL)
    checkEquals(colSds(x, na.rm=na.rm),  colSds(r, na.rm=na.rm),  msg("colSds"),  tolerance=.TOL)

    # -- mean+var combined --
    mv_r <- scRowMeanVar(x, na.rm=na.rm)
    checkEquals(mv_r[,1L], rowMeans2(r, na.rm=na.rm), msg("scRowMeanVar-mean"), tolerance=.TOL)
    checkEquals(mv_r[,2L], rowVars(r, na.rm=na.rm),   msg("scRowMeanVar-var"),  tolerance=.TOL)
    mv_c <- scColMeanVar(x, na.rm=na.rm)
    checkEquals(mv_c[,1L], colMeans2(r, na.rm=na.rm), msg("scColMeanVar-mean"), tolerance=.TOL)
    checkEquals(mv_c[,2L], colVars(r, na.rm=na.rm),   msg("scColMeanVar-var"),  tolerance=.TOL)

    # -- min / max / range --
    checkEquals(rowMins(x, na.rm=na.rm), rowMins(r, na.rm=na.rm), msg("rowMins"), tolerance=.TOL)
    checkEquals(colMins(x, na.rm=na.rm), colMins(r, na.rm=na.rm), msg("colMins"), tolerance=.TOL)
    checkEquals(rowMaxs(x, na.rm=na.rm), rowMaxs(r, na.rm=na.rm), msg("rowMaxs"), tolerance=.TOL)
    checkEquals(colMaxs(x, na.rm=na.rm), colMaxs(r, na.rm=na.rm), msg("colMaxs"), tolerance=.TOL)
    checkEquals(rowRanges(x, na.rm=na.rm), rowRanges(r, na.rm=na.rm), msg("rowRanges"), tolerance=.TOL)
    checkEquals(colRanges(x, na.rm=na.rm), colRanges(r, na.rm=na.rm), msg("colRanges"), tolerance=.TOL)

    # -- products (fully-qualified reference: the MatrixGenerics 'matrix'
    #    method leaks matrixStats' 'method' arg) --
    checkEquals(rowProds(x, na.rm=na.rm), matrixStats::rowProds(r, na.rm=na.rm), msg("rowProds"), tolerance=.TOL)
    checkEquals(colProds(x, na.rm=na.rm), matrixStats::colProds(r, na.rm=na.rm), msg("colProds"), tolerance=.TOL)

    # -- weighted means (regression for the c_colWMeans dispatch bug) --
    set.seed(100L)
    wc <- runif(ncol(x), 0.1, 2)   # length == ncol, for rowWeighted*
    wr <- runif(nrow(x), 0.1, 2)   # length == nrow, for colWeighted*
    checkEquals(rowWeightedMeans(x, wc, na.rm=na.rm), rowWeightedMeans(r, wc, na.rm=na.rm),
        msg("rowWeightedMeans"), tolerance=.TOL)
    checkEquals(colWeightedMeans(x, wr, na.rm=na.rm), colWeightedMeans(r, wr, na.rm=na.rm),
        msg("colWeightedMeans"), tolerance=.TOL)
    checkEquals(rowWeightedVars(x, wc, na.rm=na.rm), rowWeightedVars(r, wc, na.rm=na.rm),
        msg("rowWeightedVars"), tolerance=.TOL)
    checkEquals(colWeightedVars(x, wr, na.rm=na.rm), colWeightedVars(r, wr, na.rm=na.rm),
        msg("colWeightedVars"), tolerance=.TOL)

    # -- nnzero (na.counted=TRUE => NA counts as a non-zero entry) --
    checkEquals(as.numeric(row_nnzero(x, na.counted=TRUE)),
        rowSums(r != 0 | is.na(r)), msg("row_nnzero"), tolerance=.TOL)
    checkEquals(as.numeric(col_nnzero(x, na.counted=TRUE)),
        colSums(r != 0 | is.na(r)), msg("col_nnzero"), tolerance=.TOL)

    invisible(NULL)
}


# ---- example.gds : 1000 x 850 sparse counts, both orientations ----

test_stats_example <- function()
{
    fn <- system.file("extdata", "example.gds", package="SCArray")
    sce <- scExperiment(fn)
    m <- assays(sce)$counts
    checkTrue(is(m, "SC_GDSMatrix"), "counts is SC_GDSMatrix")

    r <- as.matrix(m)                   # dense reference (identical values)

    # type-1 (native, column-oriented seed)
    checkEquals(SCArray:::x_type(m), 1L, "x_type(counts)==1")
    .check_all_stats(m, r, na.rm=FALSE, tag="example/type1")

    # type-2 (transposed seed) -- same data, opposite dispatch path
    tm <- t(m)
    checkEquals(SCArray:::x_type(tm), 2L, "x_type(t(counts))==2")
    .check_all_stats(tm, t(r), na.rm=FALSE, tag="example/type2")
}


# ---- small matrix with NA values: exercises the na.rm kernels ----

test_stats_with_NA <- function()
{
    set.seed(1000L)
    nr <- 60L; nc <- 45L
    r <- matrix(0, nr, nc)
    idx <- sample.int(nr*nc, nr*nc*0.4)   # ~40% nonzero
    r[idx] <- round(runif(length(idx), 1, 50))
    r[sample.int(nr*nc, 30L)] <- NA       # scatter some NAs
    rownames(r) <- paste0("g", seq_len(nr))
    colnames(r) <- paste0("c", seq_len(nc))

    outfn <- tempfile(fileext=".gds")
    on.exit(unlink(outfn, force=TRUE))
    scConvGDS(r, outfn, verbose=FALSE)
    gf <- scOpen(outfn)
    on.exit(scClose(gf), add=TRUE, after=FALSE)
    m <- as(scArray(gf, "counts"), "SC_GDSMatrix")

    .check_all_stats(m, r, na.rm=TRUE,  tag="NA/type1/na.rm")
    .check_all_stats(m, r, na.rm=FALSE, tag="NA/type1")

    # logSumExps on safe-magnitude data (no overflow)
    rs <- r / 10
    rs[is.na(rs)] <- 0
    ms <- scReplaceNA(m, 0L) / 10
    checkEquals(unname(rowLogSumExps(ms)), unname(rowLogSumExps(rs)), "rowLogSumExps", tolerance=.TOL)
    checkEquals(unname(colLogSumExps(ms)), unname(colLogSumExps(rs)), "colLogSumExps", tolerance=.TOL)
}


# ---- useNames contract: default keeps names, matching DelayedMatrixStats ----

test_useNames <- function()
{
    fn <- system.file("extdata", "example.gds", package="SCArray")
    m <- assays(scExperiment(fn))$counts
    rn <- rownames(m); cn <- colnames(m)

    # default (useNames=TRUE): names are kept
    checkEquals(names(rowSums2(m)),  rn, "rowSums2 default names")
    checkEquals(names(colSums2(m)),  cn, "colSums2 default names")
    checkEquals(names(rowVars(m)),   rn, "rowVars default names")
    checkEquals(names(colMeans2(m)), cn, "colMeans2 default names")
    checkEquals(rownames(rowRanges(m)), rn, "rowRanges default rownames")

    # useNames=FALSE: names are dropped
    checkTrue(is.null(names(rowSums2(m, useNames=FALSE))), "rowSums2 useNames=FALSE")
    checkTrue(is.null(names(colVars(m, useNames=FALSE))),  "colVars useNames=FALSE")

    # transposed (type-2) path keeps names too
    tm <- t(m)
    checkEquals(names(rowSums2(tm)), cn, "rowSums2(t) names")
    checkEquals(names(colSums2(tm)), rn, "colSums2(t) names")
}


# ---- logSumExps numerical stability: must not overflow on large values ----

test_logSumExps_overflow <- function()
{
    set.seed(321L)
    nr <- 40L; nc <- 25L
    r <- matrix(runif(nr*nc, 500, 800), nr, nc)   # exp() would be Inf
    outfn <- tempfile(fileext=".gds")
    on.exit(unlink(outfn, force=TRUE))
    scConvGDS(r, outfn, save.sp=FALSE, verbose=FALSE)
    gf <- scOpen(outfn)
    on.exit(scClose(gf), add=TRUE, after=FALSE)
    m <- as(scArray(gf, "counts"), "SC_GDSMatrix")

    rl <- rowLogSumExps(m); cl <- colLogSumExps(m)
    checkTrue(all(is.finite(rl)), "rowLogSumExps finite (no overflow)")
    checkTrue(all(is.finite(cl)), "colLogSumExps finite (no overflow)")
    checkEquals(unname(rl), unname(matrixStats::rowLogSumExps(r)),
        "rowLogSumExps vs matrixStats", tolerance=.TOL)
    checkEquals(unname(cl), unname(matrixStats::colLogSumExps(r)),
        "colLogSumExps vs matrixStats", tolerance=.TOL)
}


# ---- rowCollapse / colCollapse (de-globalised C state) ----

test_collapse <- function()
{
    fn <- system.file("extdata", "example.gds", package="SCArray")
    m <- assays(scExperiment(fn))$counts
    r <- as.matrix(m); dimnames(r) <- NULL; dimnames(m) <- NULL
    set.seed(42L)
    ri <- sample.int(ncol(m), nrow(m), replace=TRUE)   # len == nrow
    ci <- sample.int(nrow(m), ncol(m), replace=TRUE)    # len == ncol

    checkEquals(unname(rowCollapse(m, ri)),
        unname(matrixStats::rowCollapse(r, ri)), "rowCollapse type1", tolerance=.TOL)
    checkEquals(unname(colCollapse(m, ci)),
        unname(matrixStats::colCollapse(r, ci)), "colCollapse type1", tolerance=.TOL)

    # transposed (type-2) path
    tm <- t(m); tr <- t(r)
    ri2 <- sample.int(ncol(tm), nrow(tm), replace=TRUE)
    ci2 <- sample.int(nrow(tm), ncol(tm), replace=TRUE)
    checkEquals(unname(rowCollapse(tm, ri2)),
        unname(matrixStats::rowCollapse(tr, ri2)), "rowCollapse type2", tolerance=.TOL)
    checkEquals(unname(colCollapse(tm, ci2)),
        unname(matrixStats::colCollapse(tr, ci2)), "colCollapse type2", tolerance=.TOL)

    # force many small blocks: exercises the per-block column offset logic
    old <- getAutoBlockSize(); on.exit(setAutoBlockSize(old))
    setAutoBlockSize(2e5)
    checkEquals(unname(rowCollapse(m, ri)),
        unname(matrixStats::rowCollapse(r, ri)), "rowCollapse multiblock", tolerance=.TOL)
    checkEquals(unname(colCollapse(m, ci)),
        unname(matrixStats::colCollapse(r, ci)), "colCollapse multiblock", tolerance=.TOL)
}
