Large-scale single-cell omics data manipulation using GDS files
====

![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)


## Features

SCArray provides large-scale single-cell omics data manipulation using Genomic Data Structure (GDS) files. It combines dense and sparse matrices stored in GDS files and the Bioconductor infrastructure framework (SingleCellExperiment and DelayedArray) to provide out-of-memory data storage and large-scale manipulation using the R programming language.


## Bioconductor

v1.12.0 ([http://bioconductor.org/packages/SCArray/](http://bioconductor.org/packages/SCArray/))

Package News: [NEWS](./NEWS)


## Package Maintainer

[Xiuwen Zheng](xiuwen.zheng@abbvie.com)


## Installation

* Requires R (≥ v3.5.0), [gdsfmt](http://www.bioconductor.org/packages/gdsfmt) (≥ v1.36.0)

* Bioconductor repository
```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("SCArray")
```


## Examples

```R
suppressPackageStartupMessages(library(SCArray))

# the GDS file for SingleCellExperiment
fn <- system.file("extdata", "example.gds", package="SCArray")
sce <- scExperiment(fn)

sce
## class: SingleCellExperiment
## dim: 1000 850
## metadata(0):
## assays(1): counts
## rownames(1000): MRPL20 GNB1 ... RPS4Y1 CD24
## rowData names(0):
## colnames(850): 1772122_301_C02 1772122_180_E05 ... 1772122_180_B06 1772122_180_D09
## colData names(3): Cell_ID Cell_type Timepoint
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):

counts(sce)
## <1000 x 850> sparse matrix of class SC_GDSMatrix and type "double":
##              1772122_301_C02 1772122_180_E05 1772122_300_H02 ... 1772122_180_B06
##       MRPL20               3               2               3   .               0
##         GNB1              11               6              15   .               0
##        RPL22               3               5               7   .               6
##        PARK7               1               7               3   .               2
##         ENO1               8              19              20   .               7
##          ...               .               .               .   .               .
##         SSR4               0               6               3   .               5
##        RPL10              11               4               8   .               1
## SLC25A6_loc1               4               5               5   .               3
##       RPS4Y1               0               5               0   .               2
##         CD24              18               3               7   .               0
```


## Optimized implementation of GDS-based DelayedMatrix

- ✓ = Implemented in **SCArray** for `SC_GDSMatrix`
- ⦿ = Implemented in [**DelayedMatrixStats**](http://bioconductor.org/packages/DelayedMatrixStats/) or [**MatrixGenerics**](http://bioconductor.org/packages/MatrixGenerics/)

| Method                 | Optimized | Notes                    |
|:-----------------------|:---------:|:-------------------------|
| `colAlls()`            | ⦿         |                  |
| `colAnyNAs()`          | ✓         | Check if any NA in a column  |
| `colAnys()`            | ⦿         |                  |
| `colAvgsPerRowSet()`   | ✓         | Summary statistic for equally sized subsets of rows  |
| `colCollapse()`        | ✓         | Extract one cell from each column  |
| `colCounts()`          | ⦿         |                  |
| `colCummaxs()`         | ⦿         |                  |
| `colCummins()`         | ⦿         |                  |
| `colCumprods()`        | ⦿         |                  |
| `colCumsums()`         | ⦿         |                  |
| `colDiffs()`           | ✓         | Difference between each element of a column  |
| `colIQRDiffs()`        | ⦿         |                  |
| `colIQRs()`            | ⦿         |                  |
| `colLogSumExps()`      | ✓         | Log of the sum of exponentials for each column  |
| `colMadDiffs()`        | ⦿         |                  |
| `colMads()`            | ⦿         |                  |
| `colMaxs()`            | ✓         | Maximum for each column  |
| `colMeans()`           | ✓         | Mean for each column  |
| `colMeans2()`          | ✓         | Mean for each column  |
| `colMedians()`         | ⦿         |                  |
| `colMins()`            | ✓         | Minimum for each column  |
| `colOrderStats()`      | ⦿         |                  |
| `colProds()`           | ✓         | Product for each column  |
| `colQuantiles()`       | ⦿         |                  |
| `colRanges()`          | ✓         | Minimum and maximum for each column  |
| `colRanks()`           | ⦿         |                  |
| `colSdDiffs()`         | ✓         | Standard deviation of the difference between each element of a column  |
| `colSds()`             | ✓         | Standard deviation for each column  |
| `colsum()`             | ✓         | Row sums across columns  |
| `colSums()`            | ✓         | Sum for each column  |
| `colSums2()`           | ✓         | Sum for each column  |
| `colTabulates()`       | ⦿         |                  |
| `colVarDiffs()`        | ✓         | Variance of the difference between each element of a column  |
| `colVars()`            | ✓         | Variance for each column  |
| `colWeightedMads()`    | ⦿         |                  |
| `colWeightedMeans()`   | ✓         | Weighted mean for each column  |
| `colWeightedMedians()` | ⦿         |                  |
| `colWeightedSds()`     | ✓         | Weighted standard deviation for each column  |
| `colWeightedVars()`    | ✓         | Weighted variance for each column  |
| `rowAlls()`            | ⦿         |                  |
| `rowAnyNAs()`          | ✓         | Check if any NA in a row  |
| `rowAnys()`            | ⦿         |                  |
| `rowAvgsPerColSet()`   | ✓         | Summary statistic for equally sized subsets of columns  |
| `rowCollapse()`        | ✓         | Extract one cell from each row  |
| `rowCounts()`          | ⦿         |                  |
| `rowCummaxs()`         | ⦿         |                  |
| `rowCummins()`         | ⦿         |                  |
| `rowCumprods()`        | ⦿         |                  |
| `rowCumsums()`         | ⦿         |                  |
| `rowDiffs()`           | ✓         | Difference between each element of a row  |
| `rowIQRDiffs()`        | ⦿         |                  |
| `rowIQRs()`            | ⦿         |                  |
| `rowLogSumExps()`      | ✓         | Log of the sum of exponentials for each row  |
| `rowMadDiffs()`        | ⦿         |                  |
| `rowMads()`            | ⦿         |                  |
| `rowMaxs()`            | ✓         | Maximum for each row  |
| `rowMeans()`           | ✓         | Mean for each row  |
| `rowMeans2()`          | ✓         | Mean for each row  |
| `rowMedians()`         | ⦿         |                  |
| `rowMins()`            | ✓         | Minimum for each row  |
| `rowOrderStats()`      | ⦿         |                  |
| `rowProds()`           | ✓         | Product for each row  |
| `rowQuantiles()`       | ✓         |                  |
| `rowRanges()`          | ✓         | Minimum and maximum for each column  |
| `rowRanks()`           | ⦿         |                  |
| `rowSdDiffs()`         | ✓         | Standard deviation of the difference between each element of a row  |
| `rowSds()`             | ✓         | Standard deviation for each row  |
| `rowsum()`             | ✓         | Column sums across rows  |
| `rowSums()`            | ✓         | Sum for each row  |
| `rowSums2()`           | ✓         | Sum for each row  |
| `rowTabulates()`       | ⦿         |                  |
| `rowVarDiffs()`        | ✓         | Variance of the difference between each element of a row  |
| `rowVars()`            | ✓         | Variance for each row  |
| `rowWeightedMads()`    | ⦿         |                  |
| `rowWeightedMeans()`   | ✓         | Weighted mean for each row  |
| `rowWeightedMedians()` | ⦿         |                  |
| `rowWeightedSds()`     | ✓         | Weighted standard deviation for each row  |
| `rowWeightedVars()`    | ✓         | Weighted variance for each row |
