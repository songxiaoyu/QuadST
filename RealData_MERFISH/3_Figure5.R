library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(Seurat)
library(rsq)
setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/SpatialTranscriptomics/Paper_QuadST_Revision")
load("Data/2023_merfish_frontalcortex/processed/MERFISH_scran_sce.RData")
load("Results/merfish_quadst_res.RData")


{ title_size <- 12
  legend_size <- 10
  text_size <- 11
  number_size <- 3
  axis_text_size <- 10
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
}

# Fig5a. UMAP ###################  
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


p_a1 <- ggplot(dat, aes(x=UMAP1, y=UMAP2, color = cellClass)) + 
  geom_point()+ 
  scale_color_manual(values=cbPalette, name="") + 
  theme_classic() + 
  coord_fixed()+ 
  theme(axis.title.x = element_text(size=text_size), 
        axis.title.y = element_text(size=text_size),
        axis.text.x = element_text(size=axis_text_size), 
        axis.text.y = element_text(size=axis_text_size))
p_a1

p_a2 <- ggplot(dat, aes(x=UMAP1, y=UMAP2, color = region)) + 
  geom_point() + scale_color_manual(values=brewer.pal(n = 3, name = "Purples")[1:3], name="")+
  coord_fixed()+ 
  theme_classic() + 
  theme(axis.title.x = element_text(size=text_size), 
        axis.title.y = element_text(size=text_size), 
        axis.text.x = element_text(size=axis_text_size), 
        axis.text.y = element_text(size=axis_text_size),
        legend.text=element_text(size=legend_size))# +
# guides(color = guide_legend(nrow = 2))
p_a2


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

p_b <- ggplot(dat2, aes(x=name, y=value)) + 
  geom_bar(stat="identity", width=0.7) + facet_grid(~UMAP) + 
  theme_classic() + 
  coord_fixed(ratio=6)+ 
  lims(y=c(0,1))+
  theme(axis.title.x = element_text(size=text_size), 
        axis.title.y = element_text(size=text_size), 
        axis.text.x = element_text(size=axis_text_size, angle = 60, hjust = 1), 
        axis.text.y = element_text(size=axis_text_size)) + 
  geom_text(aes(label=value), vjust=-0.2, color="blue", size=number_size) + 
  scale_x_discrete(labels=c("Cell Type", "Layer"))+ 
  ylab("Partial R-squared") + xlab("")
p_b

# Fig5c.  bubble chart ###################  

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
# levels(anchor1)=c("Astrocyte","Endothelial","Excitatory neuron","Inhibitory neuron",
#                  "Microglia","Oligodendrocyte","Oligodendrocyte\n(Precursor)","Pericyte")
neighbor1<-factor(merfish_df$neighbor)
# levels(neighbor1)=c("Astrocyte","Endothelial","Excitatory neuron","Inhibitory neuron",
#                    "Microglia","Oligodendrocyte","Oligodendrocyte\n(Precursor)","Pericyte")

merfish_df1=merfish_df
merfish_df1$anchor<-anchor1
merfish_df1$neighbor<-neighbor1

## Plot Figure
# Text sizes
{
  title_size <- 12
  legend_size <- 10
  text_size <- 11
  axis_text_size <- 10}
#seqfish_df$sig_gene_count > 0
p_c <- ggplot(merfish_df1[, ], aes(anchor,neighbor))+ geom_point(aes(size = ifelse(sig_gene_count > 0, sig_gene_count, NA)))+ scale_size_continuous(breaks=c(1, 5, 10,60, 120), limits=c(0, 300)) +
  theme_minimal() + theme(axis.title.x = element_text(size=text_size), axis.title.y = element_text(size=text_size), axis.text.x = element_text(angle = 45, hjust = 1, size=axis_text_size), 
                          axis.text.y = element_text(size=axis_text_size), plot.title = element_text(size = title_size),  legend.text=element_text(size=legend_size),legend.position = "bottom") +
  labs(size = "count") +xlab("Anchor") + ylab("Neighbor") 
p_c


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


p_d= ggplot(dat=DSres, aes(x =  neighbor, y= DS,
                              color = neighbor, alpha= ICG2)) + 
  geom_point(position = position_jitter(width = 0.1, height = 0))+
  scale_color_manual(values=cbPalette, name="Neighbor")+
  facet_grid(~ anchor, switch = "x") + 
  #coord_fixed(ratio=1)+
  scale_alpha_manual(values = c(1, 0.01))+
  theme_minimal() + 
  theme(
    strip.text.x = element_text(angle = 45),   # Customize strip text
    strip.background = element_blank(),    # Style strip background
    strip.placement = "outside",                            # Place strips outside the plot area
    axis.text.x = element_blank()
  )+ 
  labs(y="Directional Association Score", 
       x="Anchor")+
  guides(
    alpha = guide_legend(
      title = "Identification",
    ))

p_d



# Output Figure5 ===============================
library(cowplot)
p_a=p_a1 + p_a2
upper=plot_grid(p_a, p_b, labels = c("A", "B"), rel_widths = c(3.5,1), ncol = 2)
upper

lower=plot_grid(p_c,  p_d, labels = c( "C",  "D"), 
                rel_widths = c(1,1.5), ncol = 2)
lower
plot_grid(upper,  lower,  ncol=1)


ggsave(filename = "Figures/Figure5.pdf",  # File name
       width = 11,                   # Width in inches
       height = 8,                  # Height in inches
       units = "in"                # Units: "in" (inches), "cm" (centimeters), or "mm" (millimeters)
)
