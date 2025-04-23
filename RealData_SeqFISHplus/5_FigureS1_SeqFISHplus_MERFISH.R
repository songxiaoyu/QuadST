# FigureS1 cell-cell interaction distance for seqFISH+ and MERFISH    ------------
library(tidyverse)
library(ggplot2)
library(dplyr)
setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/SpatialTranscriptomics/Paper_QuadST_Revision")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# FigureS1_A cell-cell interaction distance for seqFISH+ and MERFISH    ------------
load("Results/seqfish_quadst_res.RData")
seq_df<-NULL
for(i in names(seqfish_res_list)){
  res_i<-seqfish_res_list[[i]]
  anchor_i<- strsplit(i, "-")[[1]][1]
  neighbor_i<- strsplit(i, "-")[[1]][2]
  
  seq_df$celltypepair<-c(seq_df$celltypepair, i)
  if(res_i[["summary.table"]][["sig_gene_count"]]==0){
    seq_df$distance<-c(seq_df$distance, NA)
  }else{
    seq_df$distance<-c(seq_df$distance, res_i$dist.cutoff[[1]])
  }
  
  seq_df$anchor<-c(seq_df$anchor, anchor_i)
  seq_df$neighbor<-c(seq_df$neighbor, neighbor_i)
  seq_df$group<-c(seq_df$groups, "seqFISH+")
}
seq_df<-as.data.frame(seq_df)
seq_df<-seq_df%>% group_by(anchor)%>%mutate(position = rank(-distance))

factor=(15/(102.6))
seq_df$distance_upadate<-seq_df$distance*factor

p1 <- ggplot(seq_df, aes(x = anchor, y = distance_upadate, fill = neighbor, group = position)) + 
  geom_bar(stat="Identity", position="dodge")+ 
  coord_flip() + 
  scale_fill_manual(values = cbPalette[c(2:7)], name="Neighbor") + 
  xlab("Anchor") + 
  ylab(expression(paste("Interaction distance (1", mu, "m)", sep=""))) + 
  ggtitle("SeqFISH+")+ 
  theme(axis.title = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.line = element_blank()) +
  ylim(0, 75)+
  theme_minimal()
p1

# FigureS1_B cell-cell interaction distance for MERFISH    ------------
load("Results/merfish_quadst_res.RData")

mer_df<-NULL
for(i in names(merfish_res_list)){
  res_i<-merfish_res_list[[i]]
  anchor_i<- strsplit(i, "-")[[1]][1]
  neighbor_i<- strsplit(i, "-")[[1]][2]
  
  mer_df$celltypepair=c(mer_df$celltypepair, i)
  if(res_i[["summary.table"]][["sig_gene_count"]]==0){
    mer_df$distance<-c(mer_df$distance, NA)
  }else{
    mer_df$distance<-c(mer_df$distance, res_i$dist.cutoff[[1]])
  }
  
  mer_df$anchor<-c(mer_df$anchor, anchor_i)
  mer_df$neighbor<-c(mer_df$neighbor, neighbor_i)
  mer_df$group<-c(mer_df$groups, "MERFISH")
}

mer_df<-as.data.frame(mer_df)
mer_df<-mer_df%>% group_by(anchor)%>%mutate(position = rank(-distance))
anchor1<-factor(mer_df$anchor)
levels(anchor1)=c("Astrocyte","Endothelial","Excitatory neuron","Inhibitory neuron",
                  "Microglia","Oligodendrocyte","Oligodendrocyte\n(Precursor)","Pericyte")
neighbor1<-factor(mer_df$neighbor)
levels(neighbor1)=c("Astrocyte","Endothelial","Excitatory neuron","Inhibitory neuron",
                    "Microglia","Oligodendrocyte","Oligodendrocyte\n(Precursor)","Pericyte")

mer_df$anchor<-anchor1
mer_df$neighbor<-neighbor1




p2 <- ggplot(mer_df, aes(x = anchor, y = distance, 
                         fill = neighbor, group = position)) + 
  geom_bar(stat="Identity", position="dodge") + 
  coord_flip() + 
  scale_fill_manual(values = cbPalette[c(2:9)], name="Neighbor") + 
  xlab("Anchor") + 
  ylab(expression(paste("Interaction distance (1", mu, "m)", sep=""))) + 
  ggtitle("MERFISH")+ 
  ylim(0, 75)+
  theme(axis.title = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank()) + theme_minimal()
p2

#  Output Figure3 ===============================
library(patchwork)
p= p1/ p2
p
ggsave(filename = "Figures/FigureS1.pdf",  # File name
       width = 5,                   # Width in inches
       height = 7,                  # Height in inches
       units = "in"                # Units: "in" (inches), "cm" (centimeters), or "mm" (millimeters)
)


