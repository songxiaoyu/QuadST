# Preprocess ------------------------
rm(list=ls())
library(data.table)
library(tidyverse)
library(Giotto)
library(SingleCellExperiment)
library(scran)
library(biomaRt)
library(gridExtra)
library(QuadST)
library(RANN)
library(doParallel)
library(foreach)
library(dplyr)
library(openxlsx)
setwd("QuadST")

## Read data
dat <- readRDS("Data/2023_merfish_frontalcortex/raw_data/BrainAgingSpatialAtlas_MERFISH.rds")
sce <- Seurat::as.SingleCellExperiment(dat)
# Subset 4wk donor #1 ##########################################################
# 
keepAge <- "4wk"
keepDonor <- c("MsBrainAgingSpatialDonor_1")
#c("MsBrainAgingSpatialDonor_2", "MsBrainAgingSpatialDonor_3", "MsBrainAgingSpatialDonor_4")
keepCellTypes <- c("Astro", "Endo", "ExN", "InN", "Micro", "OPC", "Olig", "Peri")
keepRegions <- c("cortical layer II/III", "cortical layer VI", "cortical layer V")

cellLocInfo <- c("center_x", "center_y", "min_x", "max_x", "min_y", "max_y", "donor_id", "age", "slice", "tissue", "fov")
cellTypeInfo <- c("cell_type", "cell_type_annot")

sce$tissue <- as.character(sce$tissue)
sce.sub <- sce[, sce$age == keepAge & sce$donor_id %in% keepDonor & sce$tissue %in% keepRegions & sce$cell_type_annot %in% keepCellTypes]
metaData.sub <- colData(sce.sub) %>% as_tibble(.)
cellAnnot <- metaData.sub %>% dplyr::select(c(cellLocInfo, cellTypeInfo))
countsData <- assay(sce.sub, "counts")
identical(colnames(countsData), rownames(colData(sce.sub)))

## Convert gene names (Ensemble to Symbol)
mouse <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", mirror = "www")
mouseAnnot <- getBM(attributes=c("ensembl_gene_id", "mgi_symbol"), filters= "ensembl_gene_id", values = rownames(countsData), mart=mouse)
mouseAnnot.ord <- mouseAnnot[match(rownames(countsData), mouseAnnot$ensembl_gene_id),]
identical(mouseAnnot.ord[,"ensembl_gene_id"], rownames(countsData))

## Rename counts data with gene symbol
rownames(countsData) <- mouseAnnot.ord[,"mgi_symbol"]

## Assign cellIDs in cellLocations, cellTypes, and counts
cellIDs <- paste0("cell", 1:dim(sce.sub)[2])
cellIDs.df <- cbind.data.frame(oldID=colnames(countsData), newID=cellIDs)
#identical(cellIDs.df[,"oldID"], colnames(countsData))
cellAnnot <- cellAnnot %>% mutate(id=cellIDs)
colnames(countsData) <- cellIDs

# Adjust  cell-specific bias factor ##########################################################
## Calculate & adjust for cell-specific bias factor
cells <- cellAnnot %>% mutate(area = abs(as.numeric(max_x) - as.numeric(min_x))*abs(as.numeric(max_y) - as.numeric(min_y))) %>% dplyr::select(id, center_x, center_y, cell_type_annot, tissue, fov, area); setnames(cells, c("id", "center_x", "center_y", "cell_type_annot", "tissue", "fov"), c("cellID", "x", "y", "cellClass", "region", "fov"))
counts <- cbind.data.frame("genes"=rownames(countsData), as.data.frame(countsData))
cells <- cells %>% group_by(cellClass)
cellTokeep <- pull(cells, cellID)
identical(colnames(counts), c("genes", cellTokeep))

## 1 Create SingleCellExperiment object ##################################
countsM <- as.data.frame(counts)
rownames(countsM) <- countsM[,1]
countsM <- countsM[,-c(1)] %>% as.matrix(.)
cellClasses <- cells %>% dplyr::select(cellID, cellClass, region, fov, area,x,y)%>%
  dplyr::mutate(x = as.numeric(as.character(x)),y = as.numeric(as.character(y)))
sce <- SingleCellExperiment(list(counts=countsM), colData=DataFrame(cellClasses), rowData=DataFrame(geneID=rownames(countsM)))

## 2 Calculate & adjust for cell-specific bias factor ##########################
sce <- computeSumFactors(sce)
sce <- logNormCounts(sce)
counts <- assay(sce, "counts")
adjusted.counts <- t(counts)/sizeFactors(sce)
adjusted.counts.t <- t(adjusted.counts)

instructions <- createGiottoInstructions(
  save_plot = FALSE,
  show_plot = TRUE,
  return_plot = FALSE,
  save_dir = 'homedir',
  python_path = NULL
)
## Create giotto object 
cell_locations1=cellClasses[,c("cellID","fov", "x", "y")]
cell_locations1 <- cell_locations1[match(colnames(adjusted.counts.t), cell_locations1$cellID), ]
MERFISH_scran <- createGiottoObject(expression =adjusted.counts.t,
                                        spatial_locs = cell_locations1,
                                        instructions = instructions)

## Add additional annotation 
cell_meta<-cellClasses[,c("cellID","cellClass","region","x", "y")]
cell_meta<-cell_meta[match(colnames(adjusted.counts.t), cell_meta$cellID), ]
cell_meta$cellClass <- droplevels(as.factor(cell_meta$cellClass))
MERFISH_scran <- addCellMetadata(MERFISH_scran,
                                     new_metadata = cell_meta,
                                     by_column = TRUE,
                                     column_cell_ID = "cellID")
# load("merfish_obj.RData")
# NCEM  ------------------------
MERFISH_scran_ncem = MERFISH_scran
cell_metadata <- pDataDT(MERFISH_scran_ncem)
uni_cell_types=levels(cell_metadata$cellClass)
cell_metadata_list <- split(cell_metadata, cell_metadata$cellClass)
res_ncem_list<-list()
for(dist_i in c(59,69,79)){
  dist.cutoff=dist_i
  expr=getExpression(gobject = MERFISH_scran_ncem,value="raw", output = "matrix")
  ## Filter high expressed genes across all cell types.
  # Filter top 75th quantile of gene expression
  # Filter gene expression with number of non zeros (at least 5 samples per gene)
  expr=as.matrix(expr)
  xm <- apply(expr, 1, mean)
  xmq <- quantile(xm, 0.25)
  gene_to_keep1 <- names(xm)[xm > xmq]
  xm2 <- apply(expr, 1, function(f) sum(f!=0)>=5)
  gene_to_keep2 <- names(xm2)[xm2]
  expr_sub<-expr[intersect(gene_to_keep1, gene_to_keep2),]
  
  # Transform cell-specific bias adjusted counts to normally distributed values.
  expr_norm <- transform_count_to_normal(expr_sub)
  
  g=nrow(expr_norm)
  n=ncol(expr_norm)
  x_l<- fastDummies::dummy_cols(cell_metadata$cellClass, remove_selected_columns=T) %>% as.matrix()
  
  compare_dist_cutoff <- function(anchor, neighbor, dist.cutoff) {
    df_anchor = cell_metadata_list[[anchor]]
    df_neighbor = cell_metadata_list[[neighbor]]
    if(anchor==neighbor){
      min_dist=sapply(1:nrow(df_anchor),function(i){
        cell <- df_anchor[i, c("x", "y")]
        dist <- sqrt((df_neighbor$x - cell[["x"]])^2 + (df_neighbor$y - cell[["y"]])^2)
        return(sort(dist)[2])
      })
    }else{
      min_dist=sapply(1:nrow(df_anchor),function(i){
        cell <- df_anchor[i, c("x", "y")]
        dist <- sqrt((df_neighbor$x - cell[["x"]])^2 + (df_neighbor$y - cell[["y"]])^2)
        return(min(dist))
      })
    }

    anchor_within_ids=df_anchor[which(min_dist<dist.cutoff),"cell_ID"]
    return(anchor_within_ids)
  }
  
  n_cores <- parallel::detectCores() - 1
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  results_list <- foreach(anchor_j = uni_cell_types, .combine = c) %dopar% {
    temp_list <- list()
    for (neighbor_i in uni_cell_types) {
      anchor_within_ids <- compare_dist_cutoff(anchor_j, neighbor_i, dist.cutoff)
      #anchor_within_ids <- anchor_within_ids[!is.na(anchor_within_ids)]
      key <- paste(anchor_j, neighbor_i, sep = "_")
      temp_list[[key]] <- anchor_within_ids
    }
    temp_list
  }
  stopCluster(cl)
  
  x_s <- data.frame(matrix(0, nrow = length(cell_metadata$cell_ID), ncol = length(uni_cell_types)))
  rownames(x_s) <- cell_metadata$cell_ID
  colnames(x_s) <- uni_cell_types
  for(anchor_j in uni_cell_types) {
    for (neighbor_i in uni_cell_types) {
      res_i<-results_list[[paste0(anchor_j,"_",neighbor_i)]]
      anchor_within_ids <- res_i
      x_s[res_i,neighbor_i]<- 1
    }
  }
  x_s=as.matrix(x_s)
  
  x_ts=sapply(1:n, function(f) matrix(outer(x_l[f,], x_s[f,]))) %>% t(.)
  tmp_n<- sub("^\\.data_", "", colnames(x_l))
  idx1=rep(tmp_n, each=8)
  idx2=rep(tmp_n, times=8)
  colnames(x_ts)=paste0("from_",idx1, "_to_", idx2)
  x_c = fastDummies::dummy_cols(cell_metadata$region, remove_selected_columns=T) %>% as.matrix()
  
  x_d = cbind(x_ts,x_l,x_c)
  y<-t(expr_norm)
  
  cl <- makeCluster(detectCores() - 1) 
  registerDoParallel(cl)  
  res_ncem_raw <- foreach(f = 1:g, .combine = 'cbind') %dopar% {
    pvals <- summary(lm(y[, f] ~ x_d - 1))$coefficients
    pvals_from <- pvals[grep("^x_dfrom", rownames(pvals)), 4]
    return(pvals_from)
  }
  stopCluster(cl)
  
  res_ncem_raw1<-t(res_ncem_raw)
  res_ncem<-data.frame(anchor=rep(idx2,each=280),neighbor=rep(idx1,each=280),
                       gene=rownames(expr_norm),pvalue=c(res_ncem_raw1))
  res_ncem<-res_ncem%>%group_by(anchor,neighbor)%>%mutate(pvalue_adj = p.adjust(pvalue,method="BH"))%>%mutate(ICG=ifelse(pvalue_adj<=0.1,1,0))
  res_ncem_list[[as.character(dist_i)]]<-res_ncem
}

# Giotto ------------------------
# Number of Exn cells is too large to directly compute distance from spatstat.geom::pairdist,
# I wrote a function separately for Exn to find minimum neighboring distance for each cell.
# Interaction Changed Features 
## Filter high expressed genes within each cell type.
# https://github.com/drieslab/Giotto/blob/a170bbeaa2ac4d4226138b8d0ac0a8cf95f4599b/R/spatial_interaction.R#L953
t.test.perGene <- function(x, i){
  if (length(x[i]) <= 1 | length(x[-i]) <= 1) {
    pvalue <- NaN 
  } else {
    pvalue <- t.test(x=x[i], y = x[-i], alternative = c("two.sided"))$p.value
  }
  return(pvalue)
}
reg.test.perGene <- function(x, i,cov=NULL){
  if (length(x[i]) <= 1 | length(x[-i]) <= 1) {
    pvalue <- NaN 
  } else {
    Indicator<-rep(0,length(x))
    Indicator[i]<-1
    if(is.null(cov)){
      pvalue<-summary(lm(x~1+Indicator))$coefficients["Indicator",4]
    }else{
      pvalue<-summary(lm(x~1+Indicator+cov))$coefficients["Indicator",4]
    }
  }
  return(pvalue)
}
compare_dist <- function(anchor, neighbor) {
  df_anchor = cell_metadata_list[[anchor]]
  df_neighbor = cell_metadata_list[[neighbor]]
  min_dist=sapply(1:nrow(df_anchor),function(i){
    cell <- df_anchor[i, c("x", "y")]
    dist <- sqrt((df_neighbor$x - cell[["x"]])^2 + (df_neighbor$y - cell[["y"]])^2)
    return(sort(dist)[2])
  })
  return(min_dist)
}
MERFISH_scran_go = MERFISH_scran
cell_metadata <- pDataDT(MERFISH_scran_go)
cell_types = unique(cell_metadata$cellClass)

expr_adj<-getExpression(gobject = MERFISH_scran_go,value="raw", output = "matrix")
res_go_list<-list()
for(dist_i in c(59,69,79)){
  res_go_raw<-NULL
  for(anchor_i in cell_types){
    for(neighbor_i in cell_types){
      anchor_idx = cell_metadata[cell_metadata$cellClass==anchor_i,cell_ID]
      anchor_cov = as.factor(cell_metadata[cell_metadata$cellClass==anchor_i,region])
      # Filter top 75th quantile of gene expression
      # Filter gene expression with number of non zeros (at least 5 samples per gene)
      expr_anchor=as.matrix(expr_adj[,c(anchor_idx)])
      xm <- apply(expr_anchor, 1, mean)
      xmq <- quantile(xm, 0.25)
      gene_to_keep1 <- names(xm)[xm > xmq]
      xm2 <- apply(expr_anchor, 1, function(f) sum(f!=0)>=5)
      gene_to_keep2 <- names(xm2)[xm2]
      expr_anchor_sel<-expr_anchor[intersect(gene_to_keep1, gene_to_keep2),]
      # Transform cell-specific bias adjusted counts to normally distributed values.
      expr_norm <- transform_count_to_normal(expr_anchor_sel)
      
      # For each anchor-neighbor cell pair, use the median distance among 
      #   the nearest neighbors as interaction cutoff.
      if(anchor_i==neighbor_i){
        if(anchor_i!="ExN"){
          meta_sel = cell_metadata[cell_metadata$cellClass==anchor_i,]
          meta_sel <- meta_sel %>% arrange(match(cell_ID, colnames(expr_norm)))
          
          d1=spatstat.geom::pairdist(meta_sel$x, meta_sel$y)
          diag(d1)=Inf
          d1min=apply(d1, 1, min)
          cutoff=dist_i
          Interact_cells=which(d1min<cutoff)
        }else{
          meta_sel = cell_metadata[cell_metadata$cellClass==anchor_i,]
          meta_sel <- meta_sel %>% arrange(match(cell_ID, colnames(expr_norm)))
          
          min_dist_i<-compare_dist(anchor_i, neighbor_i)
          cutoff=dist_i
          Interact_cells=which(min_dist_i<cutoff)
        }
      }else
      {
        meta_sel = cell_metadata[cell_metadata$cellClass==anchor_i,]
        meta_sel <- meta_sel %>% arrange(match(cell_ID, colnames(expr_norm)))
        
        meta_sel1 <- cell_metadata[cell_metadata$cellClass==neighbor_i,]
        dist_mat=spatstat.geom::crossdist.default(meta_sel$x, meta_sel$y, meta_sel1$x, meta_sel1$y)
        d1min=apply(dist_mat, 1, min)
        cutoff=dist_i
        Interact_cells=which(d1min<cutoff)
      }
      pvalue_i <- apply(expr_norm, 1, function(gene) reg.test.perGene(gene, i=Interact_cells,cov=anchor_cov) )
      pvalue_i_adj<- p.adjust(pvalue_i, method="BH")
      
      res_go_raw<-rbind(res_go_raw, data.frame(anchor = rep(anchor_i,length(pvalue_i)), neighbor = rep(neighbor_i,length(pvalue_i)), 
                                               gene=rownames(expr_norm),pvalue=pvalue_i,pvalue_adj=pvalue_i_adj,
                                               nr_select=length(Interact_cells),nr_other=length(anchor_idx)-length(Interact_cells)))
      cat(paste0(anchor_i,"-",neighbor_i),length(pvalue_i),"\n")
    }
  }
  res_go1<-res_go_raw%>%mutate(ICG=ifelse(pvalue_adj<=0.1,1,0),nr_diff=nr_select-nr_other)
  res_go<-res_go_raw%>%dplyr::select(anchor,neighbor,gene,pvalue,pvalue_adj)%>%mutate(ICG=ifelse(pvalue_adj<=0.1,1,0))
  
  res_go_list[[as.character(dist_i)]]<-res_go
}
# Giotto DT ------------------------
MERFISH_scran_go_DT = MERFISH_scran
MERFISH_scran_go_DT <- createSpatialNetwork(gobject = MERFISH_scran_go_DT,
                                                minimum_k = 2, 
                                                maximum_distance_delaunay = 10000,
                                                name="Delaunay_network")
annot_spatnetwork <- annotateSpatialNetwork(MERFISH_scran_go_DT, spatial_network_name = "Delaunay_network",cluster_column = "cellClass")
cell_metadata <- pDataDT(MERFISH_scran_go_DT)
cell_types = unique(cell_metadata$cellClass)
expr_adj<-getExpression(gobject = MERFISH_scran_go_DT,value="raw", output = "matrix")
res_go_raw<-NULL
for(anchor_i in cell_types){
  for(neighbor_i in cell_types){
    from_to_int<-paste0(neighbor_i,"-",anchor_i)
    to_from_int<-paste0(anchor_i,"-",neighbor_i)
    annot_sub=annot_spatnetwork[annot_spatnetwork$from_to%in%c(from_to_int,to_from_int),]
    anchor_cov = as.factor(cell_metadata[cell_metadata$cellClass==anchor_i,region])
    if(nrow(annot_sub)==0){
      pvalue_i=NA
      res_go_raw<-rbind(res_go_raw, data.frame(anchor = rep(anchor_i,nrow(expr_norm)), neighbor = rep(neighbor_i,nrow(expr_norm)), 
                                               gene=rownames(expr_norm),pvalue=NA,pvalue_adj=NA,nr_select=NA,nr_other=NA))
      cat(from_to_int,0,"\n")
    }else{
      anchor_idx = cell_metadata[cell_metadata$cellClass==anchor_i,cell_ID]
      # Filter top 75th quantile of gene expression
      # Filter gene expression with number of non zeros (at least 5 samples per gene)
      expr_anchor=as.matrix(expr_adj[,c(anchor_idx)])
      xm <- apply(expr_anchor, 1, mean)
      xmq <- quantile(xm, 0.25)
      gene_to_keep1 <- names(xm)[xm > xmq]
      xm2 <- apply(expr_anchor, 1, function(f) sum(f!=0)>=5)
      gene_to_keep2 <- names(xm2)[xm2]
      expr_anchor_sel<-expr_anchor[intersect(gene_to_keep1, gene_to_keep2),]
      # Transform cell-specific bias adjusted counts to normally distributed values.
      expr_norm <- transform_count_to_normal(expr_anchor_sel)
      if(anchor_i==neighbor_i){ 
        Interact_cells=match(unique(c(annot_sub$from,annot_sub$to)),colnames(expr_norm))
      }else{
        to1<-annot_sub[annot_sub$to_cell_type==anchor_i, to]
        from1<-annot_sub[annot_sub$from_cell_type==anchor_i, from]
        Interact_cells=match(unique(c(to1,from1)),colnames(expr_norm))
      }
      pvalue_i <- apply(expr_norm, 1, function(gene) reg.test.perGene(gene, i=Interact_cells,cov=anchor_cov) )
      pvalue_i_adj<- p.adjust(pvalue_i, method="BH")
      
      res_go_raw<-rbind(res_go_raw, data.frame(anchor = rep(anchor_i,length(pvalue_i)), neighbor = rep(neighbor_i,length(pvalue_i)), 
                                               gene=rownames(expr_norm),pvalue=pvalue_i,pvalue_adj=pvalue_i_adj,
                                               nr_select=length(Interact_cells),nr_other=length(anchor_idx)-length(Interact_cells)))
      cat(from_to_int,length(pvalue_i),"\n")
    }
  }
}

res_go_DT<-res_go_raw%>%dplyr::select(anchor,neighbor,gene,pvalue,pvalue_adj)%>%mutate(ICG=ifelse(pvalue_adj<=0.1,1,0))

# LR ------------------------
# Linear Regression
MERFISH_scran_lr = MERFISH_scran
cell_metadata <- pDataDT(MERFISH_scran_lr)
cell_types = levels(cell_metadata$cellClass)
expr_adj<-getExpression(gobject = MERFISH_scran_lr,value="raw", output = "matrix")
res_lr_raw<-NULL
for(anchor_i in cell_types){
  for(neighbor_i in cell_types){
    anchor_idx = cell_metadata[cell_metadata$cellClass==anchor_i,cell_ID]
    anchor_cov = as.factor(cell_metadata[cell_metadata$cellClass==anchor_i,region])
    #anchor_cov = as.factor(cell_metadata[cell_metadata$cell_types==anchor_i,FOV])
    # Filter top 75th quantile of gene expression
    # Filter gene expression with number of non zeros (at least 5 samples per gene)
    expr_anchor=as.matrix(expr_adj[,c(anchor_idx)])
    xm <- apply(expr_anchor, 1, mean)
    xmq <- quantile(xm, 0.25)
    gene_to_keep1 <- names(xm)[xm > xmq]
    xm2 <- apply(expr_anchor, 1, function(f) sum(f!=0)>=5)
    gene_to_keep2 <- names(xm2)[xm2]
    expr_anchor_sel<-expr_anchor[intersect(gene_to_keep1, gene_to_keep2),]
    # Transform cell-specific bias adjusted counts to normally distributed values.
    expr_norm <- transform_count_to_normal(expr_anchor_sel)
    
    # For each anchor-neighbor cell pair, use the median distance among 
    #   the nearest neighbors as interaction cutoff.
    if(anchor_i==neighbor_i){
      if(anchor_i!="ExN"){
        meta_sel = cell_metadata[cell_metadata$cellClass==anchor_i,]
        meta_sel <- meta_sel %>% arrange(match(cell_ID, colnames(expr_norm)))
        
        d1=spatstat.geom::pairdist(meta_sel$x, meta_sel$y)
        diag(d1)=Inf
        y=apply(d1, 1, min)
      }else{
        meta_sel = cell_metadata[cell_metadata$cellClass==anchor_i,]
        meta_sel <- meta_sel %>% arrange(match(cell_ID, colnames(expr_norm)))
        
        min_dist_i<-compare_dist(anchor_i, neighbor_i)
        y=min_dist_i
      }
    }else{
      meta_sel = cell_metadata[cell_metadata$cellClass==anchor_i,]
      meta_sel <- meta_sel %>% arrange(match(cell_ID, colnames(expr_norm)))
      meta_sel1 <- cell_metadata[cell_metadata$cellClass==neighbor_i,]
      dist_mat=spatstat.geom::crossdist.default(meta_sel$x, meta_sel$y, meta_sel1$x, meta_sel1$y)
      y=apply(dist_mat, 1, min)
    }
    pvalue_i <- sapply(1:nrow(expr_norm), function(gene) summary(lm(y ~ expr_norm[gene,anchor_idx]+anchor_cov-1))$coefficient[1,4])
    pvalue_i_adj<- p.adjust(pvalue_i, method="BH")
    
    res_lr_raw<-rbind(res_lr_raw, data.frame(anchor = rep(anchor_i,length(pvalue_i)), neighbor = rep(neighbor_i,length(pvalue_i)), 
                                             gene=rownames(expr_norm),pvalue=pvalue_i,pvalue_adj=pvalue_i_adj))
    cat(paste0(anchor_i,"-",neighbor_i),length(pvalue_i),"\n")
  }
}
res_lr<-res_lr_raw%>%mutate(ICG=ifelse(pvalue_adj<=0.1,1,0))


# save data ------------------------
merfish_ncem_list=list()
for(dist_i in c(59,69,79)){
  res_i=res_ncem_list[[as.character(dist_i)]]
  merfish_ncem_list[[as.character(dist_i)]]<-res_i%>%
    group_by(anchor,neighbor)%>%summarise(sig_gene_count=sum(ICG))%>%
    arrange(anchor, neighbor)
}
merfish_go_list=list()
for(dist_i in c(59,69,79)){
  res_i=res_go_list[[as.character(dist_i)]]
  merfish_go_list[[as.character(dist_i)]]<-res_i%>%
    group_by(anchor,neighbor)%>%summarise(sig_gene_count=sum(ICG))%>%
    arrange(anchor, neighbor)
}
merfish_go_DT=res_go_DT%>%group_by(anchor,neighbor)%>%summarise(sig_gene_count=sum(ICG))%>%arrange(anchor, neighbor)
merfish_lr=res_lr%>%dplyr::group_by(anchor,neighbor)%>%dplyr::summarise(sig_gene_count=sum(ICG))%>%dplyr::arrange(anchor, neighbor)

res_quadst=read.xlsx("SupplementaryTable2.xlsx", 
                     sheet = 2, colNames = FALSE)
colnames(res_quadst) <- as.character(res_quadst[2, ])
res_quadst <- res_quadst[-c(1,2), ]

merfish_quadst=res_quadst%>%mutate(ICG = as.numeric(ICG))%>%group_by(anchor,neighbor)%>%summarise(sig_gene_count=sum(ICG))%>%arrange(anchor, neighbor)

final_res_merfish<-cbind(rep("MERFISH",64),merfish_quadst$anchor,merfish_quadst$neighbor,
                         merfish_quadst$sig_gene_count,
                         merfish_go_DT$sig_gene_count, 
                 merfish_go_list[["59"]]$sig_gene_count,
                 merfish_go_list[["69"]]$sig_gene_count,
                 merfish_go_list[["79"]]$sig_gene_count,
                 merfish_ncem_list[["59"]]$sig_gene_count,
                 merfish_ncem_list[["69"]]$sig_gene_count,
                 merfish_ncem_list[["79"]]$sig_gene_count,
                 merfish_lr$sig_gene_count)
colnames(final_res_merfish)<-c("Dataset","Anchor","Neighbor",
                               "QuadST",
                               "Giotto_Network",
                               "Giotto (59um)","Giotto (69um)","Giotto (79um)",
                               "NCEM (59um)","NCEM (69um)","NCEM (79um)",
                               "Linear Regression")

cellClass<-table(cell_metadata$cellClass)
cellClass_num <- setNames(as.numeric(cellClass), names(cellClass))

wb_merfish <- createWorkbook()
addWorksheet(wb_merfish, "Summary")
final_res_merfish <- as.data.frame(final_res_merfish)
final_res_merfish1 <-final_res_merfish %>% mutate(across(4:ncol(.), as.integer))
writeData(wb_merfish, sheet="Summary", x=final_res_merfish1, rowNames = FALSE)

addWorksheet(wb_merfish, "Linear Regression")
writeData(wb_merfish, sheet="Linear Regression", x=res_lr, rowNames = FALSE)

addWorksheet(wb_merfish, "Giotto_Network")
writeData(wb_merfish, sheet="Giotto_Network", x=res_go_DT, rowNames = FALSE)
for(dist_i in c(59,69,79)){
  res_i=res_ncem_list[[as.character(dist_i)]]
  addWorksheet(wb_merfish, paste0("NCEM_",dist_i,"um"))
  writeData(wb_merfish, sheet=paste0("NCEM_",dist_i,"um"), x=res_i, rowNames = FALSE)
}
for(dist_i in c(59,69,79)){
  res_i=res_go_list[[as.character(dist_i)]]
  addWorksheet(wb_merfish, paste0("Giotto_",dist_i,"um"))
  writeData(wb_merfish, sheet=paste0("Giotto_",dist_i,"um"), x=res_i, rowNames = FALSE)
}
saveWorkbook(wb_merfish, "Tables/Supplementary_Table_S4.xlsx", overwrite = TRUE)





