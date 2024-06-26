---
output: github_document
---
# `QuadST: A Powerful and Robust Approach for Identifying Cell–Cell Interaction-Changed Genes on Spatially Resolved Transcriptomics`

A full description of the method can be found in our [paper](https://doi.org/10.1101/2023.12.04.570019)

## Installation

## Example of use

Below demonstrates the QuadST analysis pipeline for a single cell type pair. Multiple cell type pairs can be analyzed in parallel.


```r
library(QuadST)
library(tidyverse)
```

```
## ── Attaching packages ───────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.3.1 ──
```

```
## ✔ ggplot2 3.4.0     ✔ purrr   0.3.4
## ✔ tibble  3.1.7     ✔ dplyr   1.0.9
## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
## ✔ readr   2.1.2     ✔ forcats 0.5.1
```

```
## Warning: package 'tidyr' was built under R version 4.1.2
```

```
## Warning: package 'readr' was built under R version 4.1.2
```

```
## ── Conflicts ──────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
```

```r
library(SingleCellExperiment)
```

```
## Loading required package: SummarizedExperiment
```

```
## Loading required package: MatrixGenerics
```

```
## Loading required package: matrixStats
```

```
## 
## Attaching package: 'matrixStats'
```

```
## The following object is masked from 'package:dplyr':
## 
##     count
```

```
## 
## Attaching package: 'MatrixGenerics'
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
##     colWeightedMeans, colWeightedMedians, colWeightedSds,
##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
##     rowWeightedSds, rowWeightedVars
```

```
## Loading required package: GenomicRanges
```

```
## Loading required package: stats4
```

```
## Loading required package: BiocGenerics
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:dplyr':
## 
##     combine, intersect, setdiff, union
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
##     union, unique, unsplit, which.max, which.min
```

```
## Loading required package: S4Vectors
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following objects are masked from 'package:dplyr':
## 
##     first, rename
```

```
## The following object is masked from 'package:tidyr':
## 
##     expand
```

```
## The following objects are masked from 'package:base':
## 
##     expand.grid, I, unname
```

```
## Loading required package: IRanges
```

```
## 
## Attaching package: 'IRanges'
```

```
## The following objects are masked from 'package:dplyr':
## 
##     collapse, desc, slice
```

```
## The following object is masked from 'package:purrr':
## 
##     reduce
```

```
## Loading required package: GenomeInfoDb
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## 
## Attaching package: 'Biobase'
```

```
## The following object is masked from 'package:MatrixGenerics':
## 
##     rowMedians
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     anyMissing, rowMedians
```

```r
library(spatstat)
```

```
## Loading required package: spatstat.data
```

```
## Loading required package: spatstat.geom
```

```
## spatstat.geom 3.0-5
```

```
## 
## Attaching package: 'spatstat.geom'
```

```
## The following object is masked from 'package:SummarizedExperiment':
## 
##     shift
```

```
## The following object is masked from 'package:GenomicRanges':
## 
##     shift
```

```
## The following objects are masked from 'package:IRanges':
## 
##     reflect, shift
```

```
## Loading required package: spatstat.random
```

```
## spatstat.random 3.0-1
```

```
## Loading required package: spatstat.core
```

```
## Loading required package: nlme
```

```
## Warning: package 'nlme' was built under R version 4.1.2
```

```
## 
## Attaching package: 'nlme'
```

```
## The following object is masked from 'package:IRanges':
## 
##     collapse
```

```
## The following object is masked from 'package:dplyr':
## 
##     collapse
```

```
## Loading required package: rpart
```

```
## Warning: package 'rpart' was built under R version 4.1.2
```

```
## spatstat.core 2.4-4
```

```
## Loading required package: spatstat.linnet
```

```
## spatstat.linnet 2.3-2
```

```
## 
## spatstat 2.3-4       (nickname: 'Watch this space') 
## For an introduction to spatstat, type 'beginner'
```

```r
library(quantreg)
```

```
## Loading required package: SparseM
```

```
## 
## Attaching package: 'SparseM'
```

```
## The following object is masked from 'package:base':
## 
##     backsolve
```

```r
library(QRank)
library(ACAT)
```

### QuadST analysis pipeline
Step 1: Create an anchor-neighbor integrated matrix.

```r
data("seqFISHplus_scran_sce")
cell_id = "cellID"
cell_coord1 = "x"
cell_coord2 = "y"
cell_type = "cellClass"
anchor = "Excitatory neuron"
neighbor = "Excitatory neuron"
covariate = "FOV"

sce_an <- create_cellpair_matrix(seqFISHplus_scran_sce, cell_id, cell_coord1, cell_coord2, cell_type, anchor, neighbor, cov=covariate)
sce_an
```

```
## class: SingleCellExperiment 
## dim: 10000 325 
## metadata(0):
## assays(3): counts logcounts adjusted.counts
## rownames(10000): 1700022a21rik 1700025g04rik ... Opn1sw Pramef12
## rowData names(1): geneID
## colnames(325): cell_468 cell_469 ... cell_318 cell_304
## colData names(9): cellID cellClass ... anchor neighbor
## reducedDimNames(0):
## mainExpName: NULL
## altExpNames(0):
```

Step 2: Model and test distance-expression association

```r
# A. Determine the number of quantile levels, e.g., to ensure that there are at least 5 samples in each quantile level.
anchor_cell_count <- length(colData(sce_an)[, cell_id])
dist_taus <- create_quantile_levels(min_sample_per_quantile = 5, cell_count = anchor_cell_count, max_default = 49)

# B. Subset genes to analyze, e.g., top 25% highly expressed genes.
xm <- apply(t(assay(sce_an, "adjusted.counts")), 2, mean)
xmq <- quantile(xm, 0.75)
gene_to_keep <- names(xm)[xm > xmq]
sce_an_sub <- sce_an[gene_to_keep,]

# C. Transform cell-specific bias adjusted counts to normally distributed values.
assay(sce_an_sub, "normalized_counts") <- transform_count_to_normal(assay(sce_an_sub, "adjusted.counts"))

# D. Provide the colunm names of sce_an_sub that correspond to distance, expression values, and covariates to be used for analysis.
dist <- "distance"
expr <- "normalized_counts"
covariate <- "FOV"
QRpvalue <- test_QuadST_model(sce_an_sub, dist, expr, cov = covariate, tau = dist_taus)
str(QRpvalue)
```

```
##  num [1:2500, 1:49] 0.0983 0.5507 0.5376 0.1763 0.0205 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:2500] "Aatf" "Abat" "Abhd2" "Abhd6" ...
##   ..$ : chr [1:49] "0.02" "0.04" "0.06" "0.08" ...
```

Step 3: Identify ICGs.

```r
res <- identify_ICGs(sce_an_sub, QRpvalue, dist, expr, cov = covariate, tau = dist_taus, p_thres = 0.05, fdr = 0.1, ABconst = 0.1)

# A. Check ICGs
str(res$ICGs)
```

```
##  chr [1:303] "Adap1" "Afg3l2" "Ago2" "Aldh2" "Ap2a1" "Apba1" "Arl4d" ...
```

```r
# B. Interaction quantile and distance
str(res$q_int)
```

```
##  Named int 6
##  - attr(*, "names")= chr "0.12"
```

```r
str(res$dist_int)
```

```
##  Named num 81.4
##  - attr(*, "names")= chr "12%"
```
