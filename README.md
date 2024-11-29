
# `QuadST: A Powerful and Robust Approach for Identifying Cell–Cell Interaction-Changed Genes on Spatially Resolved Transcriptomics`

A full description of the method can be found at XXX

## Installation

With R, users can install the QuadST package directly from GitHub with
devtools:

``` r
if (!require("devtools", quietly = TRUE)) install.packages("devtools")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("SingleCellExperiment", quietly = TRUE)) BiocManager::install("SingleCellExperiment")
devtools::install_github("songxiaoyu/QuadST/Rpackage")
```

## Example of using QuadST analysis pipeline

``` r
library(QuadST)
```

### Load a public data set

Our example data was published by Eng, Chee-Huat Linus, et al 2019
(<https://www.nature.com/articles/s41586-019-1049-y>). You can
downloaded it directly in our folder Data_2019_seqfish_plus_SScortex or
using this link
(<https://github.com/songxiaoyu/QuadST/tree/main/Data_2019_seqfish_plus_SScortex>).

The data preprocessing step can be found in XXX.R, which uses three
inputs: (1) gene \* cell count matrix for expression, (2) cell spatial
coordinates, (3) cell features (e.g. cell type, FOV). The output is a
“SingleCellExperiment” object called seqFISHplus_scran_sce for QuadST
analysis. Here we directly load seqFISHplus_scran_sce to demonstrate
QuadST.

``` r
data("seqFISHplus_scran_sce")
```

#### Step 1: Specify data and parameters

Below demonstrates the QuadST analysis pipeline for an anchor-neighbor
single cell type pair (Excitatory neuron - Excitatory neuron). Multiple
cell type pairs can be analyzed in parallel.

``` r
data("seqFISHplus_scran_sce") # load a "SingleCellExperiment" object
cell_id = "cellID"  # variable name for cell index
cell_coord1 = "x" # variable name of spatial coordinate x
cell_coord2 = "y" # variable name of spatial coordinate x
cell_type = "cellClass" # variable name for cell type
anchor = "Excitatory neuron" # anchor cell type
neighbor = "Excitatory neuron" # neighbor cell type, which can be the same as anchor or different. 
datatype <- "normcounts" # data types used for QuadST analysis 
covariate <- "FOV" # covariates to adjust for 
k=1 # No. of nearest neighbors. 
d.limit=Inf # the limit of distance for cell pairing. 
```

#### Step 2: Create an anchor-neighbor integrated matrix.

``` r
sce_an <- create_integrated_matrix(seqFISHplus_scran_sce, cell_id, cell_coord1, 
                                   cell_coord2, cell_type, anchor, neighbor, k=k, d.limit = d.limit)
sce_an
```

#### Step 3: Model and test distance-expression association

``` r
# A. Determine the number of quantile levels, e.g., to ensure that there are at least 5 samples in each quantile level.

dist_taus <- create_quantile_levels(min_sample_per_quantile = 5, cell_count = dim(sce_an)[2], max_default = 49)

# B. Subset genes to analyze, e.g., top 25% highly expressed genes.
xm <- apply(t(assay(sce_an, "adjusted.counts")), 2, mean)
xmq <- quantile(xm, 0.75)
gene_to_keep <- names(xm)[xm > xmq]
sce_an_sub <- sce_an[gene_to_keep,]

# C. Transform cell-specific bias adjusted counts to normally distributed values.
assay(sce_an_sub, "normcounts") <- transform_count_to_normal(assay(sce_an_sub, "adjusted.counts"))

# D. Provide the colunm names of sce_an_sub that correspond to distance, expression values, and covariates to be used for analysis.

QRpvalue <- test_QuadST_model(x=sce_an_sub, datatype=datatype,cov = covariate, tau = dist_taus, parallel=T)
str(QRpvalue)
```

    ##  num [1:2500, 1:49] 0.78 0.675 0.969 0.789 0.725 ...
    ##  - attr(*, "dimnames")=List of 2
    ##   ..$ : chr [1:2500] "Aatf" "Abat" "Abhd2" "Abhd6" ...
    ##   ..$ : chr [1:49] "0.02" "0.04" "0.06" "0.08" ...

#### Step 4: Identify CCIs/ICGs.

``` r
res <- identify_ICGs(pMatrix=QRpvalue, fdr = 0.1)

# A. Check ICGs
res$summary.table
```

    ##   idx_ICG Q_taus sig_gene_count
    ## 1       6   0.88            296

``` r
res$data.table[1:5,] 
```

    ##        gene     pvalue      eFDR ICG
    ## Aatf   Aatf 0.11230273 0.3349762   0
    ## Abat   Abat 0.22049735 0.4651381   0
    ## Abhd2 Abhd2 0.55438487 0.6976465   0
    ## Abhd6 Abhd6 0.09536937 0.3114332   0
    ## Abl1   Abl1 0.66343834 0.7600874   0

``` r
# B. Interaction quantile and distance
distance=ICG_distance(x=sce_an, ICG.summary=res$summary.table, k=k) 
distance
```

    ##      88% 
    ## 81.40098
