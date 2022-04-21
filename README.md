Large-scale single-cell RNA-seq data manipulation using GDS files
====

![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)


## Features

Large-scale single-cell RNA-seq data manipulation and statistical analysis with scalable implementation of generalized mixed models and principal component analysis. The package integrates the sparse matrix in Genomic Data Structure (GDS) files and the Bioconductor infrastructure framework to provide out-of-memory data storage and manipulation using the R programming language.


## Bioconductor:

v1.4.0 ([http://bioconductor.org/packages/SCArray/](http://bioconductor.org/packages/SCArray/))

Package News: [NEWS](./NEWS)


## Package Maintainer

[Xiuwen Zheng](xiuwen.zheng@abbvie.com)


## Installation

* Requires R (≥ v3.5.0), [gdsfmt](http://www.bioconductor.org/packages/gdsfmt) (≥ v1.24.0)

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
fn <- system.file("extdata", "LaMannoBrainData.gds", package="SCArray")
sce <- scExperiment(fn)

sce
## class: SingleCellExperiment
## dim: 12000 243
## metadata(0):
## assays(1): counts
## rownames(12000): Rp1 Sox17 ... Efhd2 Fhad1
## rowData names(0):
## colnames(243): 1772072122_A04 1772072122_A05 ... 1772099011_H05 1772099012_E04
## colData names(2): CELL_ID Cell_type
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):

counts(sce)
## <12000 x 243> sparse matrix of class SC_GDSMatrix and type "double":
##         1772072122_A04 1772072122_A05 1772072122_A06 ... 1772099011_H05 1772099012_E04
##     Rp1              0              0              0   .              0              0
##   Sox17              0              0              0   .              0              0
##  Mrpl15              1              2              1   .              2              2
##  Lypla1              0              0              1   .              0              1
##   Tcea1              1              0              4   .              6              1
##     ...              .              .              .   .              .              .
##   Agmat              0              0              0   .              0              0
## Dnajc16              0              0              0   .              0              0
##   Casp9              0              0              0   .              0              0
##   Efhd2              0              0              0   .              1              1
##   Fhad1              0              0              0   .              1              0
```

