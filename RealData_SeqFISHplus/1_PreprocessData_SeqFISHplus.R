# Title ##################################################################
## Purpose: Preprocess seqFISHplus data for QuadST analysis 
## Analysis:
##
## Input: 1. cortex_svz_expression.txt
##        2. cortex_svz_centroids_coord.txt
##        3. cortex_svz_centroids_annot.txt
## Output: A list of "SingleCellExperiment" objects containing all cell types
## "seqFISHplus_scran_sce_list" for QuadST analysis
##
# Preprocess Data ###################################################################
rm(list=ls())
library(data.table)
library(tidyverse)
library(Giotto)
library(SingleCellExperiment)
library(scran)
library(gridExtra)

setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/SpatialTranscriptomics/Paper_QuadST_Revision")

## load data
SS_locations <- data.table::fread("Data/2019_seqfish_plus_SScortex/cell_locations/cortex_svz_centroids_coord.txt")
cortex_fields <- data.table::fread("Data/2019_seqfish_plus_SScortex/cell_locations/cortex_svz_centroids_annot.txt")
counts <- data.table::fread("Data/2019_seqfish_plus_SScortex/count_matrix/cortex_svz_expression.txt") 

# Create a stitch file ##########################################################
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

# Subset cortex and celltypes ##########################################################
## Subset cells and counts data from cortex and celltypes
cells <- cell_locations; setnames(cells, c("ID", "X", "Y", "cell_types", "X_final", "Y_final"), 
                                  c("cellID", "X_raw", "Y_raw", "cellClass", "x", "y"))
setnames(counts, "V1", "genes")

FOV_cortex <- c("FOV #0", "FOV #1", "FOV #2", "FOV #3", "FOV #4")
cellClass <- table(cells$cellClass)
cells <- cells %>% filter(FOV %in% FOV_cortex) %>% group_by(cellClass) %>% filter(n() >= 10)
cellTokeep <- pull(cells, cellID)
## https://cloud.r-project.org/web/packages/data.table/vignettes/datatable-intro.html
## Why "with = FALSE" is necessary for data table
counts <- counts[, c("genes", cellTokeep), with = FALSE] %>% as.data.frame()
cellClasses <- cells %>% dplyr::select(cellID, cellClass, FOV,x,y)

# Adjust  cell-specific bias factor ##########################################################
## Calculate & adjust for cell-specific bias factor 
## 1 Create SingleCellExperiment object ##################################
countsM <- as.data.frame(counts)
rownames(countsM) <- countsM[,1]
countsM <- countsM[,-c(1)] %>% as.matrix(.)

cellClasses <- cells %>% dplyr::select(cellID, cellClass, FOV,x,y)

seqFISHplus_scran_sce <- SingleCellExperiment(list(counts=countsM), 
                                              colData=DataFrame(cellClasses), 
                                              rowData=DataFrame(geneID=rownames(countsM)),
                                              metadata="2019_seqfish_plus_SScortex")

## 2 Calculate & adjust for cell-specific bias factor ################
seqFISHplus_scran_sce <- computeSumFactors(seqFISHplus_scran_sce)
seqFISHplus_scran_sce <- logNormCounts(seqFISHplus_scran_sce)
counts <- assay(seqFISHplus_scran_sce, "counts")
adjusted.counts <- t(counts)/sizeFactors(seqFISHplus_scran_sce)
assay(seqFISHplus_scran_sce, "adjusted.counts") <- t(adjusted.counts)


# Save SingleCellExperiment objects list ##############################
save(seqFISHplus_scran_sce,  file="Data/2019_seqfish_plus_SScortex/processed/seqFISHplus_scran_sce.RData")
