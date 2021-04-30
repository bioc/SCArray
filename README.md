Large-scale single-cell RNA-seq data manipulation using GDS files
====

![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)
[GNU General Public License, GPLv3](http://www.gnu.org/copyleft/gpl.html)


## Features

Large-scale single-cell RNA-seq data manipulation and statistical analysis with scalable implementation of generalized mixed models and principal component analysis. The package integrates the sparse matrix in Genomic Data Structure (GDS) files and the Bioconductor infrastructure framework to provide out-of-memory data storage and manipulation using the R programming language.


## Bioconductor:

v1.0.0 ([http://bioconductor.org/packages/SCArray/](http://bioconductor.org/packages/SCArray/))


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
