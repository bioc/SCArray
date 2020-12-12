suppressPackageStartupMessages(library(RUnit))
suppressPackageStartupMessages(library(SCArray))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(DelayedMatrixStats))


test_sce_matrix <- function()
{
	# a GDS file for SingleCellExperiment
	fn <- system.file("extdata", "LaMannoBrainData.gds", package="SCArray")

	sce <- scExperiment(fn)
	mt <- assays(sce)$counts

	cM <- unname(colMeans(mt))
	cM2 <- colMeans2(mt)
	checkEquals(cM, cM2, "counts: colmean")
	checkEquals(min(cM), 0.1414167, "counts: colmean min", tolerance=1e-6)
	checkEquals(max(cM), 1.689167, "counts: colmean max", tolerance=1e-6)

	rM <- unname(rowMeans(mt))
	rM2 <- rowMeans2(mt)
	checkEquals(rM, rM2, "counts: rowmean")
	checkEquals(min(rM), 0, "counts: rowmean min", tolerance=1e-6)
	checkEquals(max(rM), 477.1934, "counts: rowmean max", tolerance=1e-6)
}


test_sce_col_row <- function()
{
	# a GDS file for SingleCellExperiment
	fn <- system.file("extdata", "LaMannoBrainData.gds", package="SCArray")
	sce <- scExperiment(fn)

	# column data
	cD <- colData(sce)
	checkEquals(dim(cD), c(243, 2), "column data: dim")

	# row data
	rD <- rowData(sce)
	checkEquals(dim(rD), c(12000, 0), "row data: dim")
}
