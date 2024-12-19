rm(list=ls())
library(QuadST)
library(parallel)
library(tidyverse)
library(SingleCellExperiment)
setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/SpatialTranscriptomics/Paper_QuadST_Revision")

# Load Data------------------
load("Data/2019_seqfish_plus_SScortex/processed/seqFISHplus_scran_sce.RData")
# Run QuadST ------------------
{ # fixed parameter
  cell_id = "cellID"
  cell_coord1 = "x"
  cell_coord2 = "y"
  cell_type = "cellClass"
  # anchor = "Excitatory neuron"
  # neighbor = "Interneuron"
  cov <- "FOV"
  k=1
  d.limit=Inf}

RunQuadST=function(anchor, neighbor){
  # create SingleCellExperiment object 
  sce_an <- QuadST::create_integrated_matrix(seqFISHplus_scran_sce, 
                                             cell_id=cell_id, cell_coord1, cell_coord2, cell_type, 
                                             anchor=anchor, neighbor=neighbor, k, d.limit)
  # Determine the number of quantile levels, e.g., 
  # to ensure that there are at least 5 samples in each quantile level.
  dist_taus <- QuadST::create_quantile_levels(min_sample_per_quantile = 5, cell_count = dim(sce_an)[2], 
                                              max_default = 49)
  
  # Filter top 25th quantile of gene expression
  # Filter gene expression with number of non zeros (at least 5 samples per gene)
  expr=assay(sce_an, "adjusted.counts")
  xm <- apply(expr, 1, mean)
  xmq <- quantile(xm, 0.75)
  gene_to_keep1 <- names(xm)[xm > xmq]
  xm2 <- apply(expr, 1, function(f) sum(f!=0)>5)
  gene_to_keep2 <- names(xm2)[xm2]
  sce_an_sub <- sce_an[intersect(gene_to_keep1, gene_to_keep2),]
  
  # Transform cell-specific bias adjusted counts to normally distributed values.
  assay(sce_an_sub, "normcounts") <- transform_count_to_normal(assay(sce_an_sub, "adjusted.counts"))
  
  # Provide the colunm names of sce_an_sub that correspond to distance, expression values, and covariates to be used for analysis.
  
  QRpvalue <- QuadST::test_QuadST_model(x=sce_an_sub, datatype="normcounts", 
                                        cov = cov, tau = dist_taus, parallel=T)
  res <- identify_ICGs(pMatrix=QRpvalue, fdr = 0.1)
  
  #get distance
  distance<-QuadST::ICG_distance(x=sce_an, ICG.summary=res$summary.table, k=k)
  res$dist.cutoff<-distance
  res$distance<-data.frame(cellID=sce_an@colData@listData[["cellID"]],
                           w_distance=sce_an@colData@listData[["w_distance"]])
  
  # get Directional Score
  tau=res$summary.table[2] %>% as.numeric()
  y = sce_an_sub$w_distance
  x<-sce_an_sub@assays@data@listData$normcounts
  covariate=sce_an_sub$FOV

  coef=sapply(1:nrow(x), function(f) quantreg::rq(y ~ x[f,]+ covariate,tau=tau)$coefficients['x[f, ]'])
  res$data.table=data.frame(res$data.table, coef)
  return(res)
}

celltypes=unique(seqFISHplus_scran_sce$cellClass)
seqfish_res_list=NULL
for (i in celltypes) {
  for (j in celltypes){
    cell_pair=paste0(i,"-",j)  
    cat(cell_pair,"\n")
    res_i=RunQuadST(anchor=i, neighbor=j) 
    seqfish_res_list[[cell_pair]]=res_i
  }
}

save(seqfish_res_list,file="Results/seqfish_quadst_res.RData")


# output csv as supplementary table --------



res <- NULL
for(i in 1:length(seqfish_res_list)){
  name<-names(seqfish_res_list)[i]
  anchor_i<- strsplit(name, "-")[[1]][1]
  neighbor_i<- strsplit(name, "-")[[1]][2]
  sig_gene_i<-seqfish_res_list[[i]]$data.table
  res<-rbind(res, data.frame(anchor = anchor_i, neighbor = neighbor_i, sig_gene_i))
}
write.table(res, file="Result/SupplementaryTable1.csv", sep=",", row.names = F, col.names = T)

