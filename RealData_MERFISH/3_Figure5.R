# Figure5-merfish. pdf: A+B+C+D
## *---------------------------
##
## Script name: Figure5
##
## Purpose of script: spatial map(5a) + umap(5b) + partial r2(5c) +bubble chart(5d)
## + direction score(5e)
##
## *---------------------------
##
## Notes:
##   
## *---------------------------
rm(list=ls())
library(RColorBrewer)
library(ggplot2)
library(tidyr)
library(dplyr)
library(Seurat)
library(rsq)
setwd("QuadST")
load("Results/merfish_quadst_res.RData")
load("Data/2023_merfish_frontalcortex/processed/MERFISH_scran_sce.RData")

{ title_size <- 12
  legend_size <- 10
  text_size <- 11
  number_size <- 3
  axis_text_size <- 11
  cbPalette <- c("#999999", "#E69F00", "#009E73", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
}

# Fig5a. Spatial map ###################  
#p_a
dt_mer<-MERFISH_scran_sce@colData %>%as.data.frame()
dt_mer$cellClass<-factor(dt_mer$cellClass)
levels(dt_mer$cellClass)=c("Astrocyte","Endothelial","Excitatory neuron","Inhibitory neuron",
                           "Microglia","Oligodendrocyte","Oligodendrocyte\n(Precursor)","Pericyte")

## Plot Figure 3a: p_a-----------------
p_a<-ggplot(dt_mer, aes(x = x, y = y, color = cellClass)) +
  geom_point( size = 0.1) +  
  scale_color_manual(values=cbPalette, name="") +
  labs(x = "X coords",
       y = "Y coords") +
  coord_fixed()+
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size = 1))) +
  theme( legend.position = "right",        
         legend.justification = c(0, 1),    
         legend.background = element_rect(fill = NA, color = NA),
         axis.title.x = element_text(size = text_size),  
         axis.title.y = element_text(size = text_size),  
         axis.text.x = element_text(size=axis_text_size), 
         axis.text.y = element_text(size=axis_text_size),
         legend.text=element_text(size=legend_size),
         legend.title=element_text(size=text_size),
         plot.margin = margin(t = 10, r = 10, b = 10, l = 0) )
p_a



# Fig5b. UMAP ###################  
obj_seurat <- as.Seurat(MERFISH_scran_sce, counts = "counts")
obj_seurat <- NormalizeData(obj_seurat)
obj_seurat <- FindVariableFeatures(obj_seurat)
obj_seurat <- ScaleData(obj_seurat)
obj_seurat <- RunPCA(obj_seurat)
obj_seurat <- RunUMAP(obj_seurat, dims = 1:10)
DimPlot(obj_seurat, reduction = "umap", group.by = "cellClass")
dat=data.frame(obj_seurat@reductions$umap@cell.embeddings[,1:2], 
               obj_seurat@meta.data)
colnames(dat)[1:2]=c("UMAP1", "UMAP2")

dat$cellClass<-factor(dat$cellClass)
levels(dat$cellClass)=c("Astrocyte","Endothelial","Excitatory neuron","Inhibitory neuron",
                           "Microglia","Oligodendrocyte","Oligodendrocyte\n(Precursor)","Pericyte")


p_b1 <- ggplot(dat, aes(x=UMAP1, y=UMAP2, color = cellClass)) + 
  geom_point()+ 
  scale_color_manual(values=cbPalette, name="") + 
  theme_classic() + 
  coord_fixed()+ 
  theme(axis.title.x = element_text(size=text_size), 
        axis.title.y = element_text(size=text_size),
        axis.text.x = element_text(size=axis_text_size), 
        axis.text.y = element_text(size=axis_text_size),
        legend.text=element_text(size=legend_size),
        legend.position = "top",  
        legend.title = element_blank(),
        legend.key.size = unit(1, 'pt'))+
  guides(color = guide_legend(nrow = 3))#+
  #guides(color = guide_legend(nrow = 2))
p_b1

dat$region <- gsub("cortical layer ", "", dat$region)
p_b2 <- ggplot(dat, aes(x=UMAP1, y=UMAP2, color = region)) + 
  geom_point() + scale_color_manual(values=brewer.pal(n = 3, name = "Purples")[1:3], name = "Cortical Layer")+
  coord_fixed()+ 
  theme_classic() + 
  theme(axis.title.x = element_text(size=text_size), 
        axis.title.y = element_text(size=text_size), 
        axis.text.x = element_text(size=axis_text_size), 
        axis.text.y = element_text(size=axis_text_size),
        legend.text=element_text(size=legend_size),
        legend.position = "top",  
        legend.margin = margin(t = 5, r = 5, b = 0, l = 5)) 
p_b2


# Fig5b. Partial R2 ###################  


# umap_1
fit1 <- lm(UMAP1 ~  region + cellClass, data=dat)
partial.rsq1 <- rsq.partial(fit1)$partial.rsq %>% round(., digits=2)

fit2 <- lm(UMAP2 ~  region + cellClass, data=dat)
partial.rsq2 <- rsq.partial(fit2)$partial.rsq  %>% round(., digits=2)

effects=rbind(partial.rsq1, partial.rsq2) %>% as.data.frame() %>%
  mutate(UMAP=c("UMAP1", "UMAP2")) %>%
  rename(Layer=V1, CellType=V2) 


dat2=effects %>% pivot_longer(cols = 1:2)

p_c <- ggplot(dat2, aes(x=name, y=value)) + 
  geom_bar(stat="identity", width=0.7) + facet_grid(~UMAP) + 
  theme_classic() + 
  coord_fixed(ratio=6)+ 
  lims(y=c(0,1))+
  theme(axis.title.x = element_text(size=text_size), 
        axis.title.y = element_text(size=text_size), 
        axis.text.x = element_text(size=axis_text_size, angle = 55, hjust = 1), 
        axis.text.y = element_text(size=axis_text_size)) + 
  geom_text(aes(label=value), vjust=-0.2, color="blue", size=number_size) + 
  scale_x_discrete(labels=c("Cell Type", "Layer"))+ 
  ylab("Partial R-squared") + xlab("")
p_c

# Fig5d.  bubble chart ###################  

results_list <- NULL
for(i in 1:length(merfish_res_list)){
  name<-names(merfish_res_list)[i]
  anchor_i<- strsplit(name, "-")[[1]][1]
  neighbor_i<- strsplit(name, "-")[[1]][2]
  sig_gene_count_i<-merfish_res_list[[i]]$summary.table$sig_gene_count
  
  results_list[[i]] <- data.frame(anchor = anchor_i, neighbor = neighbor_i, 
                                  sig_gene_count = sig_gene_count_i)
}
merfish_df <- do.call(rbind, results_list)

anchor1<-factor(merfish_df$anchor)
levels(anchor1)=c("Astrocyte","Endothelial","Excitatory neuron","Inhibitory neuron",
                 "Microglia","Oligodendrocyte","Oligodendrocyte\n(Precursor)","Pericyte")
neighbor1<-factor(merfish_df$neighbor)
levels(neighbor1)=c("Astrocyte","Endothelial","Excitatory neuron","Inhibitory neuron",
                   "Microglia","Oligodendrocyte","Oligodendrocyte\n(Precursor)","Pericyte")

merfish_df1=merfish_df
merfish_df1$anchor<-anchor1
merfish_df1$neighbor<-neighbor1

## Plot Figure
# Text sizes
#seqfish_df$sig_gene_count > 0
p_d <- ggplot(merfish_df1[, ], aes(anchor,neighbor))+ 
  geom_point(aes(size = ifelse(sig_gene_count > 0, sig_gene_count, NA)))+ 
  scale_size_continuous(breaks=c(1, 10, 150,300), limits=c(0, 3000)) +
  theme_minimal() + 
  theme(axis.title.x = element_text(size=10), 
        axis.title.y = element_text(size=10), 
        axis.text.x = element_text(angle = 45, hjust = 1, size=10), 
        axis.text.y = element_text(size=10), 
        plot.title = element_text(size = 10),  
        legend.text=element_text(size=10),
        legend.position = "bottom") +
  labs(size = "count") +xlab("Anchor") + ylab("Neighbor") 
p_d


# Fig5d.  DS ###################  

DSres <- NULL
for(i in 1:length(merfish_res_list)){
  name<-names(merfish_res_list)[i]
  anchor_i<- strsplit(name, "-")[[1]][1]
  neighbor_i<- strsplit(name, "-")[[1]][2]
  
  t<-merfish_res_list[[i]]$data.table
  DS=sign(t$coef) * (-log10(t$pvalue))
  DSres[[i]] <- data.frame(anchor = anchor_i, neighbor = neighbor_i, 
                           gene=t$gene, ICG=t$ICG, DS = DS)
}
DSres <- do.call(rbind, DSres)
DSres$ICG2=ifelse(DSres$ICG==1, "ICG", "Not ICG")

DSres$anchor<-factor(DSres$anchor)
levels(DSres$anchor)=c("Astrocyte","Endothelial","Excitatory neuron","Inhibitory neuron",
                  "Microglia","Oligodendrocyte","Oligodendrocyte\n(Precursor)","Pericyte")
DSres$neighbor<-factor(DSres$neighbor)
levels(DSres$neighbor)=c("Astrocyte","Endothelial","Excitatory neuron","Inhibitory neuron",
                    "Microglia","Oligodendrocyte","Oligodendrocyte\n(Precursor)","Pericyte")


p_e= ggplot(dat=DSres, aes(x =  neighbor, y= DS,
                              color = neighbor, alpha= ICG2)) + 
  geom_point(position = position_jitter(width = 0.1, height = 0))+
  scale_color_manual(values=cbPalette, name="Neighbor")+
  facet_grid(~ anchor, switch = "x") + 
  #coord_fixed(ratio=1)+
  scale_alpha_manual(values = c(1, 0.01))+
  theme_minimal() + 
  theme(
    strip.text.x = element_text(angle = 20,size=10),   # Customize strip text
    strip.background = element_blank(),    # Style strip background
    strip.placement = "outside",                            # Place strips outside the plot area
    legend.title=element_text(size=11),
    legend.text=element_text(size=10),
    axis.title.x = element_text(size=12), 
    axis.title.y = element_text(size=12), 
    axis.text.x = element_blank()
  )+ 
  labs(y="Directional Association Score", 
       x="Anchor")+
  guides(
    alpha = guide_legend(
      title = "Identification",
      override.aes = list(alpha =c(1, 0.2))
    ))

p_e



# Output Figure5 ===============================
library(cowplot)
p_b=plot_grid(p_b1, p_b2, ncol = 2) 
upper=plot_grid(p_a, p_d, labels = c("A", "D"), rel_widths = c(1.5, 1), ncol = 2)
upper
middle=plot_grid( p_b, p_c, labels = c("B", "C"), rel_widths = c(4, 1 ), ncol = 2)
middle
lower=plot_grid(p_e, labels = c("E"),  ncol = 1)

lower
plot_grid(upper,  middle, lower,  ncol=1)

ggsave(filename = "Figures/Figure5.pdf",  # File name
       width = 11,                   # Width in inches
       height = 12,                  # Height in inches
       units = "in"                # Units: "in" (inches), "cm" (centimeters), or "mm" (millimeters)
)
