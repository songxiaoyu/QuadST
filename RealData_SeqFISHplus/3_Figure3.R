# Figure3-seqfish. pdf: scatter + umap + partial R^2
## *---------------------------
##
## Script name: Figure3 
##
## Purpose of script: spatial_scatter_plot (a) + UMAP (b) + partial_R^2 for seqFISH+ data.
##
## *---------------------------
##
## Notes:
##   
## *---------------------------
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyverse)
library(readxl)
library(WriteXLS)
library(ggrepel)
library(ggpubr)
# library(ComplexHeatmap)
library(egg)
library(rsq)
library(RColorBrewer)
library(gridExtra)
setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/SpatialTranscriptomics/Paper_QuadST_Revision")
load("Data/2019_seqfish_plus_SScortex/processed/seqFISHplus_scran_sce.RData")

## common figure setting ========================
{ title_size <- 12
  legend_size <- 10
  text_size <- 11
  number_size <- 3
  axis_text_size <- 10
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
}

# Plot Fig3a: spatial map ===============================
dt_seq<-seqFISHplus_scran_sce@colData


p_a<-ggplot(dt_seq, aes(x = x, y = y, color = cellClass)) +
  geom_point(alpha = 0.8, size = 2) +  
  scale_color_manual(values=cbPalette[2:8], name="") +
  coord_cartesian(xlim = c(min(dt_seq$x)*0.8, max(dt_seq$x)*1), 
                  ylim =c(min(dt_seq$y)*1.05, max(dt_seq$y)*1.2)) +  
  labs(x = "X coords",
       y = "Y coords") +
  coord_fixed()+
  theme_classic() +
  theme(legend.position = "right",  legend.justification = c(0, 1),
        legend.background = element_rect(fill = NA, color = NA),
        axis.title.x = element_text(size=text_size), 
        axis.title.y = element_text(size=text_size), 
        axis.text.x = element_text(size=axis_text_size), 
        axis.text.y = element_text(size=axis_text_size))
p_a
# Plot Fig3a: UMAP (b)  ===============================
seqFISHplus_seurat <- as.Seurat(seqFISHplus_scran_sce, counts = "counts")
seqFISHplus_seurat <- NormalizeData(seqFISHplus_seurat)
seqFISHplus_seurat <- FindVariableFeatures(seqFISHplus_seurat)
seqFISHplus_seurat <- ScaleData(seqFISHplus_seurat)
seqFISHplus_seurat <- RunPCA(seqFISHplus_seurat)
seqFISHplus_seurat <- RunUMAP(seqFISHplus_seurat, dims = 1:10)
DimPlot(seqFISHplus_seurat, reduction = "umap", group.by = "cellClass")
dat=data.frame(seqFISHplus_seurat@reductions$umap@cell.embeddings[,1:2], 
                                 seqFISHplus_seurat@meta.data)
colnames(dat)[1:2]=c("UMAP1", "UMAP2")


p_b1 <- ggplot(dat, aes(x=UMAP1, y=UMAP2, color = cellClass)) + 
  geom_point()+ 
  scale_color_manual(values=cbPalette[2:8], name="") + 
  theme_classic() + 
  coord_fixed()+ 
  theme(axis.title.x = element_text(size=text_size), 
        axis.title.y = element_text(size=text_size),
        axis.text.x = element_text(size=axis_text_size), 
        axis.text.y = element_text(size=axis_text_size))
p_b1

p_b2 <- ggplot(dat, aes(x=UMAP1, y=UMAP2, color = FOV)) + 
  geom_point() + scale_color_manual(values=brewer.pal(n = 6, name = "Purples")[2:6], name="")+
  coord_fixed()+ 
  theme_classic() + 
  theme(axis.title.x = element_text(size=text_size), 
        axis.title.y = element_text(size=text_size), 
        axis.text.x = element_text(size=axis_text_size), 
        axis.text.y = element_text(size=axis_text_size),
        legend.text=element_text(size=legend_size))# +
  # guides(color = guide_legend(nrow = 2))
p_b2

## Plot Figure 3c: p_c=========================


# umap_1
fit1 <- lm(UMAP1 ~  FOV + cellClass, data=dat)
partial.rsq1 <- rsq.partial(fit1)$partial.rsq %>% round(., digits=2)

fit2 <- lm(UMAP2 ~  FOV + cellClass, data=dat)
partial.rsq2 <- rsq.partial(fit2)$partial.rsq  %>% round(., digits=2)

effects=rbind(partial.rsq1, partial.rsq2) %>% as.data.frame() %>%
  mutate(UMAP=c("UMAP1", "UMAP2")) %>%
  rename(FOV=V1, CellType=V2) 
  

dat2=effects %>% pivot_longer(cols = 1:2)

p_c <- ggplot(dat2, aes(x=name, y=value)) + 
  geom_bar(stat="identity", width=0.7) + facet_grid(~UMAP) + 
  theme_classic() + 
  coord_fixed(ratio=6)+ 
  lims(y=c(0,1))+
  theme(axis.title.x = element_text(size=text_size), 
        axis.title.y = element_text(size=text_size), 
        axis.text.x = element_text(size=axis_text_size, angle = 60, hjust = 1), 
        axis.text.y = element_text(size=axis_text_size)) + 
  geom_text(aes(label=value), vjust=-0.2, color="blue", size=number_size) + 
  scale_x_discrete(labels=c("Cell Type", "FOV"))+ 
  ylab("Partial R-squared") + xlab("")
p_c

# Output Figure3 ===============================
library(cowplot)
upper=plot_grid(p_a, p_c, labels = c("A", "C"), rel_widths = c(2.4,1), ncol = 2)
lower=plot_grid(p_b1, p_b2,  labels = c( "B"), rel_widths = c(1,1), ncol = 2)
plot_grid(upper, lower,  ncol=1)


ggsave(filename = "Figures/Figure3.pdf",  # File name
       width = 11,                   # Width in inches
       height = 7,                  # Height in inches
       units = "in"                # Units: "in" (inches), "cm" (centimeters), or "mm" (millimeters)
)

