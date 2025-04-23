# TableS3-seqfish. xlsx: Summary sheet, Linear Regression, Giotto Network, NCEM(59um), NCEM(69um), NCEM(79um), Giotto(59um), Giotto(69um), Giotto(79um)
## *---------------------------
##
## Script name: TableS3
##
## Purpose of script: algorithm comparison
##
## *---------------------------
##
## Notes:
##   
## *---------------------------
rm(list=ls())
library(QuadST)
library(data.table)
library(tidyverse)
library(Giotto)
library(SingleCellExperiment)
library(scran)
library(gridExtra)
library(foreach)
library(doParallel)
library(openxlsx)
library(conflicted)
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("colnames", "base")
conflicts_prefer(base::intersect)
conflicts_prefer(base::setdiff)
conflicts_prefer(base::rownames)
#setwd("/QuadST")

## load data
SS_locations <- data.table::fread("Data/2019_seqfish_plus_SScortex/cell_locations/cortex_svz_centroids_coord.txt")
cortex_fields <- data.table::fread("Data/2019_seqfish_plus_SScortex/cell_locations/cortex_svz_centroids_annot.txt")
counts <- data.table::fread("Data/2019_seqfish_plus_SScortex/count_matrix/cortex_svz_expression.txt") 
expr_path<-"Data/2019_seqfish_plus_SScortex/count_matrix/cortex_svz_expression.txt"
## 1.1 Create a stitch file ####################################
# This dataset contains multiple field of views, which need to be stitched together

## 1 first merge location and additional metadata
SS_loc_annot <- data.table::merge.data.table(SS_locations, cortex_fields, by = 'ID')
SS_loc_annot[, ID := factor(ID, levels = paste0('cell_',1:nrow(SS_loc_annot)))]
data.table::setorder(SS_loc_annot, ID) # reindex

## 2 create file with offset information 
my_offset_file <- data.table::data.table(field = c(0, 1, 2, 3, 4, 5, 6),
                                         x_offset = c(0, 1654.97, 1750.75, 1674.35, 675.5, 2048, 675),
                                         y_offset = c(0, 0, 0, 0, -1438.02, -1438.02, 0))

## 3 create a stitch file
## stitchFieldCoordinates: giotto function
cell_locations <- stitchFieldCoordinates(location_file = SS_loc_annot,
                                         offset_file = my_offset_file,
                                         cumulate_offset_x = T,
                                         cumulate_offset_y = F,
                                         field_col = 'FOV',
                                         reverse_final_x = F,
                                         reverse_final_y = T) %>% 
  mutate(FOV=paste0("FOV #", FOV)) %>% mutate(ID = as.character(ID))

## 1.2 Subset cortex and celltypes ####################################
## Subset cells and counts data from cortex and celltypes
cells <- cell_locations; setnames(cells, c("ID", "X", "Y", "cell_types", "X_final", "Y_final"), 
                                  c("cellID", "X_raw", "Y_raw", "cellClass", "x", "y"))
setnames(counts, "V1", "genes")

FOV_cortex <- c("FOV #0", "FOV #1", "FOV #2", "FOV #3", "FOV #4")
cellClass <- table(cells$cellClass)
cells <- cells %>%as.data.frame()%>% filter(FOV %in% FOV_cortex) %>% group_by(cellClass) %>% filter(n() >= 10)
cellTokeep <- pull(cells, cellID)
## https://cloud.r-project.org/web/packages/data.table/vignettes/datatable-intro.html
## Why "with = FALSE" is necessary for data table
counts <- counts[, c("genes", cellTokeep), with = FALSE] %>% as.data.frame()
cellClasses <- cells %>% dplyr::select(cellID, cellClass, FOV,x,y)

## 1.3 Adjust  cell-specific bias factor ###########################################
instructions <- createGiottoInstructions(
  save_plot = FALSE,
  show_plot = TRUE,
  return_plot = FALSE,
  save_dir = 'homedir',
  python_path = NULL
)
countsM <- as.data.frame(counts)
rownames(countsM) <- countsM[,1]
countsM <- countsM[,-c(1)] %>% as.matrix(.)
## Create giotto object 
cell_locations1=cellClasses%>%ungroup()%>%dplyr::select(cellID, x, y) %>%dplyr::arrange(match(cellID, cellTokeep))
seqFISHplus_scran <- createGiottoObject(expression =countsM,
                                             spatial_locs = cell_locations1,
                                             instructions = instructions)
## Calculate & adjust for cell-specific bias factor 
### Add adjusted expression mat; Follow SingleCellExperiment adjust method.
load("Data/2019_seqfish_plus_SScortex/processed/seqFISHplus_scran_sce.RData")
adj_expr=assay(seqFISHplus_scran_sce,"adjusted.counts") # t(expr)/sizeFactors(seqFISHplus_scran_sce)
tmp=createExprObj(adj_expr)
seqFISHplus_scran=setExpression(seqFISHplus_scran,tmp,name="adjusted") # spat_unit="cell",feat_type="rna",

## Add additional annotation 
cortex_fields1<-cortex_fields%>%filter(ID %in% cellTokeep)%>% 
  left_join(cells[, c("cellID","x", "y")], by = c("ID"="cellID"))%>%
  arrange(match(ID, cellTokeep))
cortex_fields_update<-cortex_fields1
factor=(15/(102.6)) #0.1461988

cortex_fields_update$x=cortex_fields1$x*factor
cortex_fields_update$y=cortex_fields1$y*factor

cortex_fields_update$cell_types<-as.factor(cortex_fields_update$cell_types)
seqFISHplus_scran <- addCellMetadata(seqFISHplus_scran,
                              new_metadata = cortex_fields_update, #cortex_fields1
                              by_column = TRUE,
                              column_cell_ID = "ID")
# NCEM  ------------------------
seqFISHplus_scran_ncem = seqFISHplus_scran
cell_metadata <- pDataDT(seqFISHplus_scran_ncem)
res_ncem_list<-list()
for(dist_i in c(59,69,79)){
  dist.cutoff=dist_i
  expr=getExpression(gobject = seqFISHplus_scran_ncem,value="adjusted", output = "matrix")
  ## Filter high expressed genes across all cell types.
  # Filter top 25th quantile of gene expression
  # Filter gene expression with number of non zeros (at least 5 samples per gene)
  expr=as.matrix(expr)
  xm <- apply(expr, 1, mean)
  xmq <- quantile(xm, 0.75)
  gene_to_keep1 <- names(xm)[xm > xmq]
  xm2 <- apply(expr, 1, function(f) sum(f!=0)>=5)
  gene_to_keep2 <- names(xm2)[xm2]
  expr_sub<-expr[intersect(gene_to_keep1, gene_to_keep2),]
  
  # Transform cell-specific bias adjusted counts to normally distributed values.
  expr_norm <- transform_count_to_normal(expr_sub)
  
  g=nrow(expr_norm)
  n=ncol(expr_norm)
  y=t(expr_norm)
  x_l= fastDummies::dummy_cols(cell_metadata$cell_types, remove_selected_columns=T) %>% as.matrix()
  x_c = fastDummies::dummy_cols(cell_metadata$FOV, remove_selected_columns=T) %>% as.matrix()
  
  dist_mat=spatstat.geom::pairdist(X=cell_metadata[, 4:5])
  diag(dist_mat)<-Inf
  A_mat=1*(dist_mat<dist.cutoff)
  x_s=1*(A_mat %*% x_l>0)
  
  x_ts=sapply(1:n, function(f) matrix(outer(x_l[f,], x_s[f,]))) %>% t(.)
  tmp_n<- sub("^\\.data_", "", colnames(x_l))
  idx1=rep(tmp_n, each=6)
  idx2=rep(tmp_n, times=6)
  colnames(x_ts)=paste0("from_",idx1, "_to_", idx2)
  coef_names<-paste0("x_d",colnames(x_ts))

  y<-t(expr_norm)
  x_d = cbind(x_ts,x_l,x_c)
  
  cl <- makeCluster(detectCores() - 1) 
  registerDoParallel(cl)  
  res_ncem_raw <- foreach(f = 1:g, .combine = 'cbind') %dopar% {
    pvals=summary(lm(y[,f] ~ x_d - 1))$coefficients[ ,4]
    missing_coef=setdiff(coef_names,names(pvals))
    zero_vec <- setNames(rep("NA", length(missing_coef)), missing_coef)
    pvals1 <- c(pvals, zero_vec)
    pvals1 <- pvals1[coef_names]
    return(pvals1)
  }
  stopCluster(cl)
  res_ncem_raw1<-t(res_ncem_raw)
  res_ncem<-data.frame(anchor=rep(idx2,each=2500),neighbor=rep(idx1,each=2500),
                       gene=rownames(expr_norm),pvalue=c(res_ncem_raw1))
  res_ncem<-res_ncem%>%group_by(anchor,neighbor)%>%mutate(pvalue_adj = p.adjust(pvalue,method="BH"))%>%mutate(ICG=ifelse(pvalue_adj<=0.1,1,0))
  res_ncem_list[[as.character(dist_i)]]<-res_ncem
}

# Giotto ------------------------
# Interaction Changed Features 
## Filter high expressed genes within each cell type.
# https://github.com/drieslab/Giotto/blob/a170bbeaa2ac4d4226138b8d0ac0a8cf95f4599b/R/spatial_interaction.R#L953
t.test.perGene <- function(x, i){
  if (length(x[i]) <= 1 | length(x[-i]) <= 1) {
    pvalue <- NaN 
  } else{
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
seqFISHplus_scran_go = seqFISHplus_scran
cell_metadata <- pDataDT(seqFISHplus_scran_go)
cell_types = levels(cell_metadata$cell_types)

expr_adj<-getExpression(gobject = seqFISHplus_scran_go,value="adjusted", output = "matrix")
res_go_list<-list()
for(dist_i in c(59,69,79)){
  res_go_raw<-NULL
  for(anchor_i in cell_types){
    for(neighbor_i in cell_types){
      anchor_idx = cell_metadata[cell_metadata$cell_types==anchor_i,cell_ID]
      anchor_cov = as.factor(cell_metadata[cell_metadata$cell_types==anchor_i,FOV])
      # Filter top 25th quantile of gene expression
      # Filter gene expression with number of non zeros (at least 5 samples per gene)
      expr_anchor=as.matrix(expr_adj[,c(anchor_idx)])
      xm <- apply(expr_anchor, 1, mean)
      xmq <- quantile(xm, 0.75)
      gene_to_keep1 <- names(xm)[xm > xmq]
      xm2 <- apply(expr_anchor, 1, function(f) sum(f!=0)>=5)
      gene_to_keep2 <- names(xm2)[xm2]
      expr_anchor_sel<-expr_anchor[base::intersect(gene_to_keep1, gene_to_keep2),]
      # Transform cell-specific bias adjusted counts to normally distributed values.
      expr_norm <- transform_count_to_normal(expr_anchor_sel)
      
      # For each anchor-neighbor cell pair, use the median distance among
      #   the nearest neighbors as interaction cutoff.
      if(anchor_i==neighbor_i){
        meta_sel = cell_metadata[cell_metadata$cell_types==anchor_i,]
        meta_sel <- meta_sel %>% arrange(match(cell_ID, colnames(expr_norm)))
        
        d1=spatstat.geom::pairdist(meta_sel$x, meta_sel$y)
        diag(d1)=Inf
        d1min=apply(d1, 1, min)
        cutoff=dist_i
        Interact_cells=which(d1min<cutoff)
      }else{
        meta_sel = cells[cells$cellClass==anchor_i,]
        meta_sel <- meta_sel %>% arrange(match(cellID, colnames(expr_norm)))
        
        meta_sel1 <- cells[cells$cellClass==neighbor_i,]
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
  res_go<-res_go_raw%>%select(anchor,neighbor,gene,pvalue,pvalue_adj)%>%mutate(ICG=ifelse(pvalue_adj<=0.1,1,0))
  res_go_list[[as.character(dist_i)]]<-res_go
}
# Giotto DT ################################
seqFISHplus_scran_go_DT = seqFISHplus_scran
seqFISHplus_scran_go_DT <- createSpatialNetwork(gobject = seqFISHplus_scran_go_DT,
                                             minimum_k = 2, 
                                             maximum_distance_delaunay = 10000,
                                             name="Delaunay_network")
annot_spatnetwork <- annotateSpatialNetwork(seqFISHplus_scran_go_DT, spatial_network_name = "Delaunay_network",cluster_column = "cell_types")
cell_metadata <- pDataDT(seqFISHplus_scran_go_DT)
cell_types = levels(cell_metadata$cell_types)
expr_adj<-getExpression(gobject = seqFISHplus_scran_go_DT,value="adjusted", output = "matrix")
res_go_raw<-NULL
for(anchor_i in cell_types){
  for(neighbor_i in cell_types){
    from_to_int<-paste0(neighbor_i,"-",anchor_i)
    to_from_int<-paste0(anchor_i,"-",neighbor_i)
    annot_sub=annot_spatnetwork[annot_spatnetwork$from_to%in%c(from_to_int,to_from_int),]
    anchor_cov = as.factor(cell_metadata[cell_metadata$cell_types==anchor_i,FOV])
    if(nrow(annot_sub)==0){
      pvalue_i="NA"
      res_go_raw<-rbind(res_go_raw, data.frame(anchor = rep(anchor_i,nrow(expr_norm)), neighbor = rep(neighbor_i,nrow(expr_norm)), 
                                               gene=rownames(expr_norm),pvalue="NA",pvalue_adj="NA",nr_select="NA",nr_other="NA"))
      cat(from_to_int,0,"\n")
    }else{
      anchor_idx = cell_metadata[cell_metadata$cell_types==anchor_i,cell_ID]
      # Filter top 25th quantile of gene expression
      # Filter gene expression with number of non zeros (at least 5 samples per gene)
      expr_anchor=as.matrix(expr_adj[,c(anchor_idx)])
      xm <- apply(expr_anchor, 1, mean)
      xmq <- quantile(xm, 0.75)
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
seqFISHplus_scran_lr = seqFISHplus_scran
cell_metadata <- pDataDT(seqFISHplus_scran_lr)
cell_types = levels(cell_metadata$cell_types)
expr_adj<-getExpression(gobject = seqFISHplus_scran_lr,value="adjusted", output = "matrix")
res_lr_raw<-NULL
for(anchor_i in cell_types){
  for(neighbor_i in cell_types){
    anchor_idx = cell_metadata[cell_metadata$cell_types==anchor_i,cell_ID]
    anchor_cov = as.factor(cell_metadata[cell_metadata$cell_types==anchor_i,FOV])
    # Filter top 25th quantile of gene expression
    # Filter gene expression with number of non zeros (at least 5 samples per gene)
    expr_anchor=as.matrix(expr_adj[,c(anchor_idx)])
    xm <- apply(expr_anchor, 1, mean)
    xmq <- quantile(xm, 0.75)
    gene_to_keep1 <- names(xm)[xm > xmq]
    xm2 <- apply(expr_anchor, 1, function(f) sum(f!=0)>=5)
    gene_to_keep2 <- names(xm2)[xm2]
    expr_anchor_sel<-expr_anchor[intersect(gene_to_keep1, gene_to_keep2),]
    # Transform cell-specific bias adjusted counts to normally distributed values.
    expr_norm <- transform_count_to_normal(expr_anchor_sel)
    
    # For each anchor-neighbor cell pair, use the median distance among 
    #   the nearest neighbors as interaction cutoff.
    if(anchor_i==neighbor_i){
      meta_sel = cell_metadata[cell_metadata$cell_types==anchor_i,]
      meta_sel <- meta_sel %>% arrange(match(cell_ID, colnames(expr_norm)))
      
      d1=spatstat.geom::pairdist(meta_sel$x, meta_sel$y)
      diag(d1)=Inf
      y=apply(d1, 1, min)
    }else{
      meta_sel = cells[cells$cellClass==anchor_i,]
      meta_sel <- meta_sel %>% arrange(match(cellID, colnames(expr_norm)))
      meta_sel1 <- cells[cells$cellClass==neighbor_i,]
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
seqfish_ncem_list=list()
for(dist_i in c(59,69,79)){
  res_i=res_ncem_list[[as.character(dist_i)]]
  seqfish_ncem_list[[as.character(dist_i)]]<-res_i%>%
    group_by(anchor,neighbor)%>%summarise(sig_gene_count=sum(ICG))%>%
    arrange(anchor, neighbor)
}
seqfish_go_list=list()
for(dist_i in c(59,69,79)){
  res_i=res_go_list[[as.character(dist_i)]]
  seqfish_go_list[[as.character(dist_i)]]<-res_i%>%
    group_by(anchor,neighbor)%>%summarise(sig_gene_count=sum(ICG))%>%
    arrange(anchor, neighbor)
}
seqfish_lr=res_lr%>%group_by(anchor,neighbor)%>%summarise(sig_gene_count=sum(ICG))%>%arrange(anchor, neighbor)
seqfish_go_DT=res_go_DT%>%group_by(anchor,neighbor)%>%summarise(sig_gene_count=sum(ICG))%>%arrange(anchor, neighbor)
res_quadst=read.xlsx("SupplementaryTable1.xlsx", 
                     sheet = 2, colNames = FALSE)
colnames(res_quadst) <- as.character(res_quadst[2, ])
res_quadst <- res_quadst[-c(1,2), ]

seqfish_quadst=res_quadst%>%group_by(anchor,neighbor)%>%summarise(sig_gene_count=sum(ICG))%>%arrange(anchor, neighbor)

cellClass<-table(cells$cellClass)
cellClass_num <- setNames(as.numeric(cellClass), names(cellClass))

final_res_seqfish<-cbind(rep("SeqFish+",36),seqfish_quadst$anchor,seqfish_quadst$neighbor,
                 seqfish_quadst$sig_gene_count,
                 seqfish_go_DT$sig_gene_count, 
                 seqfish_go_list[["59"]]$sig_gene_count,
                 seqfish_go_list[["69"]]$sig_gene_count,
                 seqfish_go_list[["79"]]$sig_gene_count,
                 seqfish_ncem_list[["59"]]$sig_gene_count,
                 seqfish_ncem_list[["69"]]$sig_gene_count,
                 seqfish_ncem_list[["79"]]$sig_gene_count,
                 seqfish_lr$sig_gene_count)

colnames(final_res_seqfish)<-c("Dataset","Anchor","Neighbor",
                               "QuadST",
                               "Giotto_Network",
                               "Giotto (59um)","Giotto (69um)","Giotto (79um)",
                               "NCEM (59um)","NCEM (69um)","NCEM (79um)",
                               "Linear Regression")

wb_seqfish <- createWorkbook()
addWorksheet(wb_seqfish, "Summary")
final_res_seqfish <- as.data.frame(final_res_seqfish)
final_res_seqfish1 <-final_res_seqfish %>% mutate(across(4:ncol(.), as.integer))
writeData(wb_seqfish, sheet="Summary", x=final_res_seqfish1, rowNames = FALSE)

addWorksheet(wb_seqfish, "Linear Regression")
writeData(wb_seqfish, sheet="Linear Regression", x=res_lr, rowNames = FALSE)

addWorksheet(wb_seqfish, "Giotto_Network")
writeData(wb_seqfish, sheet="Giotto_Network", x=res_go_DT, rowNames = FALSE)
for(dist_i in c(59,69,79)){
  res_i=res_ncem_list[[as.character(dist_i)]]
  addWorksheet(wb_seqfish, paste0("NCEM_",dist_i,"um"))
  writeData(wb_seqfish, sheet=paste0("NCEM_",dist_i,"um"), x=res_i, rowNames = FALSE)
}
for(dist_i in c(59,69,79)){
  res_i=res_go_list[[as.character(dist_i)]]
  addWorksheet(wb_seqfish, paste0("Giotto_",dist_i,"um"))
  writeData(wb_seqfish, sheet=paste0("Giotto_",dist_i,"um"), x=res_i, rowNames = FALSE)
}

saveWorkbook(wb_seqfish, "Tables/Supplementary_Table_S3.xlsx", overwrite = TRUE)
