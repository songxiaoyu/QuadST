# Title ##################################################################
## Purpose: Preprocess MERFISH data for QuadST analysis 
## Analysis:
##
## Input: BrainAgingSpatialAtlas_MERFISH.rds
## Output: A "SingleCellExperiment" object MERFISH_scran_sce for QuadST analysis
## Note: BrainAgingSpatialAtlas_MERFISH.rds can be downloaed directly from
## https://cellxgene.cziscience.com/collections/31937775-0602-4e52-a799-b6acdd2bac2e



rm(list=ls())
library(data.table)
library(tidyverse)
library(Giotto)
library(SingleCellExperiment)
library(scran)
library(biomaRt)
library(gridExtra)
library(QuadST)
setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/SpatialTranscriptomics/Paper_QuadST_Revision")

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
assay(sce, "adjusted.counts") <- t(adjusted.counts)
MERFISH_scran_sce<-sce


# Save SingleCellExperiment objects list ##############################
save(MERFISH_scran_sce,file="Data/2023_merfish_frontalcortex/processed/MERFISH_scran_sce.RData")








