
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
sessionInfo() 
```

    ## R version 4.4.1 (2024-06-14)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Sonoma 14.4
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: Asia/Singapore
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] QuadST_0.1.0                SingleCellExperiment_1.26.0
    ##  [3] SummarizedExperiment_1.34.0 Biobase_2.64.0             
    ##  [5] GenomicRanges_1.56.2        GenomeInfoDb_1.40.1        
    ##  [7] IRanges_2.38.1              S4Vectors_0.42.1           
    ##  [9] BiocGenerics_0.50.0         MatrixGenerics_1.16.0      
    ## [11] matrixStats_1.4.1           BiocManager_1.30.25        
    ## [13] devtools_2.4.5              usethis_3.0.0              
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] xfun_0.49               htmlwidgets_1.6.4       remotes_2.5.0          
    ##  [4] processx_3.8.4          lattice_0.22-6          callr_3.7.6            
    ##  [7] vctrs_0.6.5             tools_4.4.1             ps_1.8.1               
    ## [10] curl_6.0.1              Matrix_1.7-0            desc_1.4.3             
    ## [13] lifecycle_1.0.4         GenomeInfoDbData_1.2.12 compiler_4.4.1         
    ## [16] MatrixModels_0.5-3      QRank_1.0               SparseM_1.84-2         
    ## [19] httpuv_1.6.15           quantreg_5.99.1         htmltools_0.5.8.1      
    ## [22] yaml_2.3.10             later_1.3.2             crayon_1.5.3           
    ## [25] urlchecker_1.0.1        MASS_7.3-61             ellipsis_0.3.2         
    ## [28] cachem_1.1.0            DelayedArray_0.30.1     sessioninfo_1.2.2      
    ## [31] abind_1.4-8             mime_0.12               digest_0.6.37          
    ## [34] purrr_1.0.2             splines_4.4.1           fastmap_1.2.0          
    ## [37] grid_4.4.1              cli_3.6.3               SparseArray_1.4.8      
    ## [40] magrittr_2.0.3          S4Arrays_1.4.1          survival_3.7-0         
    ## [43] pkgbuild_1.4.4          UCSC.utils_1.0.0        promises_1.3.0         
    ## [46] rmarkdown_2.29          XVector_0.44.0          httr_1.4.7             
    ## [49] memoise_2.0.1           shiny_1.9.1             evaluate_1.0.1         
    ## [52] knitr_1.49              miniUI_0.1.1.1          profvis_0.4.0          
    ## [55] rlang_1.1.4             Rcpp_1.0.13-1           xtable_1.8-4           
    ## [58] glue_1.8.0              pkgload_1.4.0           rstudioapi_0.17.1      
    ## [61] jsonlite_1.8.9          R6_2.5.1                fs_1.6.5               
    ## [64] zlibbioc_1.50.0

### Load a public data set

Our example data was published by Eng, Chee-Huat Linus, et al 2019
(<https://www.nature.com/articles/s41586-019-1049-y>). You can
downloaded it directly in our folder Data_2019_seqfish_plus_SScortex or
using this link
(<https://github.com/songxiaoyu/QuadST/tree/main/Data_2019_seqfish_plus_SScortex>).

The data preprocessing step can be found in
RealData_SeqFISHplus/1_PreprocessData_SeqFISHplus.R, which uses three
inputs: (1) gene \* cell count matrix for expression, (2) cell spatial
coordinates, (3) cell features (e.g. cell type, FOV). The output is a
“SingleCellExperiment” object called seqFISHplus_scran_sce for QuadST
analysis. Here we directly load seqFISHplus_scran_sce to demonstrate
QuadST.

``` r
data("seqFISHplus_scran_sce") # load a "SingleCellExperiment" object
```

#### Step 1: Specify data and parameters

Below demonstrates the QuadST analysis pipeline for an anchor-neighbor
single cell type pair (Excitatory neuron - Excitatory neuron). Multiple
cell type pairs can be analyzed in parallel.

``` r
cell_id = "cellID"  # variable name for cell index
cell_coord1 = "x" # variable name of spatial coordinate x
cell_coord2 = "y" # variable name of spatial coordinate x
cell_type = "cellClass" # variable name for cell type
anchor = "Excitatory neuron" # anchor cell type
neighbor = "Excitatory neuron" # neighbor cell type, which can be the same as anchor or different. 
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

    ## class: SingleCellExperiment 
    ## dim: 10000 325 
    ## metadata(1): ''
    ## assays(3): counts logcounts adjusted.counts
    ## rownames(10000): 1700022a21rik 1700025g04rik ... Opn1sw Pramef12
    ## rowData names(1): geneID
    ## colnames(325): cell_1 cell_10 ... cell_98 cell_99
    ## colData names(9): cellID cellClass ... anchor neighbor
    ## reducedDimNames(0):
    ## mainExpName: NULL
    ## altExpNames(0):

#### Step 3: Model and test distance-expression association

``` r
# A. Determine the number of quantile levels, e.g., to ensure that there are at 
# least 5 samples in each quantile level and a max o f49 quantile levels.

dist_taus <- create_quantile_levels(cell_count = dim(sce_an)[2], min_sample_per_quantile = 5, 
                                    max_default = 49)

# B. Subset genes to analyze, e.g., top 25% highly expressed genes.
xm <- apply(t(assay(sce_an, "adjusted.counts")), 2, mean)
xmq <- quantile(xm, 0.75)
gene_to_keep <- names(xm)[xm > xmq]
sce_an_sub <- sce_an[gene_to_keep,]

# C. Transform adjusted counts to normally distributed values with minmimal zero.
assay(sce_an_sub, "normcounts") <- transform_count_to_normal(assay(sce_an_sub, "adjusted.counts"))

# D. Provide the colunm names of sce_an_sub that correspond to distance, expression values, and covariates to be used for analysis.

QRpvalue <- test_QuadST_model(x=sce_an_sub, datatype="normcounts",cov = covariate, tau = dist_taus, 
                              parallel=T)
QRpvalue[1:3,1:3]
```

    ##             0.02      0.04      0.06
    ## Aatf  0.07556242 0.1234996 0.1127361
    ## Abat  0.20801442 0.2882366 0.3112828
    ## Abhd2 0.49058778 0.6451283 0.2081084

#### Step 4: Identify CCIs/ICGs.

``` r
res <- identify_ICGs(pMatrix=QRpvalue, fdr = 0.1)

# A. Check ICGs
res$summary.table
```

    ##   idx_ICG Q_taus sig_gene_count
    ## 1      11   0.22            326

``` r
res$data.table[1:5,] 
```

    ##        gene       pvalue      eFDR ICG
    ## Aatf   Aatf 0.2323383146 0.3761063   0
    ## Abat   Abat 0.0002044411 0.2089304   0
    ## Abhd2 Abhd2 0.7203990108 0.8399361   0
    ## Abhd6 Abhd6 0.0571624488 0.3819104   0
    ## Abl1   Abl1 0.3999580098 0.5967196   0

``` r
# B. Interaction quantile and distance
distance=ICG_distance(x=sce_an, ICG.summary=res$summary.table, k=k) 
distance
```

    ##      22% 
    ## 92.91461
