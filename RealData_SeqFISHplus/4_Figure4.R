rm(list=ls())

library(ggplot2)
library(dplyr)
library(ggrepel)
setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/SpatialTranscriptomics/Paper_QuadST_Revision")
load("Results/seqfish_quadst_res.RData")
names(seqfish_res_list)

{title_size <- 12
  legend_size <- 10
  text_size <- 11
  axis_text_size <- 10}

# Fig4a. seqfish - bubble chart ###################  


res <- NULL
for(i in 1:length(seqfish_res_list)){
  name<-names(seqfish_res_list)[i]
  anchor_i<- strsplit(name, "-")[[1]][1]
  neighbor_i<- strsplit(name, "-")[[1]][2]
  sig_gene_count_i<-seqfish_res_list[[i]]$summary.table$sig_gene_count
  res[[i]] <- data.frame(anchor = anchor_i, neighbor = neighbor_i, 
                                  sig_gene_count = sig_gene_count_i)
}
seqfish_df <- do.call(rbind, res)

# Text sizes
fig4_a <- ggplot(seqfish_df[, ], aes(anchor,neighbor))+ 
  geom_point(aes(size = ifelse(sig_gene_count > 0, sig_gene_count, NA)))+ 
  scale_size_continuous(breaks=c(1, 5, 10,60, 120), limits=c(0, 300)) +
  theme_minimal() + 
  coord_fixed(ratio=0.5)+
  theme(axis.title.x = element_text(size=text_size), 
        axis.title.y = element_text(size=text_size), 
        axis.text.x = element_text(angle = 45, hjust = 1, size=axis_text_size), 
        axis.text.y = element_text(size=axis_text_size), 
        plot.title = element_text(size = title_size),  
        legend.text=element_text(size=legend_size),
        legend.position = "bottom") +
  labs(size = "count") +xlab("Anchor") + ylab("Neighbor") 
fig4_a


# Fig4b. composition chart ###################  
TF=readxl::read_excel("Data/Database/41467_2017_BFncomms15089_MOESM3449_ESM.xlsx", sheet="allTFs")[,1] %>%as.matrix() %>%  as.character()
LR=read.table("Data/Database/mouse_lr_pair.txt", header=T)[,2:3] %>% as.matrix() %>%  as.character()

ICG=NULL
for(i in 1:length(seqfish_res_list)){
  dat=seqfish_res_list[[i]]$data.table
  ICG_i<-dat[which(dat$ICG==1), "gene"]
  ICG=c(ICG, ICG_i)
}
ICG=unique(ICG)

a=length(ICG)
# [1] 571
b=length(intersect(TF, ICG))
c=length(intersect(LR, ICG))
dat=data.frame(cat=c("Ligand & Receptor", "Transcription factor & cofactor", "Others"), 
               prop=c(b/a, c/a, (a-b-c)/a)*100)
dat$cat <- factor(dat$cat, levels = c("Ligand & Receptor", "Transcription factor & cofactor", "Others"))

fig4_b=ggplot(data=dat,aes(x = "", y = prop, fill = cat))+
  geom_bar(stat = "identity", position = "fill")+
  scale_y_continuous(labels = scales::percent_format()) +
  labs( x = NULL,
    y = "Percentage (%) ",
    fill = " "
  ) +
  coord_fixed(ratio=2.5)+
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.position = "bottom")+
  scale_fill_manual(values = c( "#00CCFF", "#FFCC00",  "#999999"))+
  guides(fill = guide_legend(nrow = 3)) 

fig4_b

# Fig4c. pathway enrichment for EN-EN ###################  
mouseGOCC_database <- GSA::GSA.read.gmt("Data/Database/m5.go.cc.v2024.1.Mm.symbols.gmt")
mouseGOCC <- mouseGOCC_database$genesets
names(mouseGOCC) <- mouseGOCC_database$geneset.names

dat<-seqfish_res_list[["Excitatory neuron-Excitatory neuron"]]$data.table
dat$logq<-(-log10(dat$eFDR))
dat$rank=rank(dat$logq)

path_p=NULL
for (i in 1:length(mouseGOCC)){
  path1=mouseGOCC[[i]]
  G_in=intersect(dat$gene, path1)
  if (length(G_in) > 5) {
      print(i)
      Q_in = dat[G_in, "rank"]
      Q_out = dat[setdiff(dat$gene, G_in), "rank"]
      test <- wilcox.test(Q_in, Q_out)
      p1=data.frame(name=names(mouseGOCC)[[i]], p=test$p.value)
      path_p=rbind(path_p, p1)
  }
}
path_p$BH=p.adjust(path_p$p, method="BH")
path_p$color="Others"
path_p$color[grep("SYNA", path_p$name)]="Neuron Signaling"
path_p[path_p$name=="GOCC_NEURON_PROJECTION",]$color="Neuron Signaling"

dat1=path_p[which(path_p$BH<0.1), ] 

transformed_strings <- sapply(dat1$name, function(x) {
  parts <- unlist(strsplit(x, "_")) # Split string by '_'
  paste(parts[-1], collapse = " ") # Remove first element and combine the rest
}); dat1$ID=transformed_strings

dat1$logp= (-log10(dat1$p))

dat1=dat1%>% arrange(logp) %>%
  mutate(ID = factor(ID, levels = ID))

fig4_c=ggplot(dat1, aes(x=ID, y= logp)) +
  geom_bar(stat = "identity")+
  coord_flip()+
  theme_classic() +
  labs(x="GO Cellular Component", y= "-log10(p-value)")+
  theme(legend.title  = element_blank())
fig4_c


# Fig4d. Directional Score ###################  

DSres <- NULL
for(i in 1:length(seqfish_res_list)){
  name<-names(seqfish_res_list)[i]
  anchor_i<- strsplit(name, "-")[[1]][1]
  neighbor_i<- strsplit(name, "-")[[1]][2]
  
  t<-seqfish_res_list[[i]]$data.table
  DS=sign(t$coef) * (-log10(t$pvalue))
  DSres[[i]] <- data.frame(anchor = anchor_i, neighbor = neighbor_i, 
                         gene=t$gene, ICG=t$ICG, DS = DS)
}
DSres <- do.call(rbind, DSres)
DSres$ICG2=ifelse(DSres$ICG==1, "ICG", "Not ICG")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

fig4_d= ggplot(dat=DSres, aes(x =  neighbor, y= DS,
                      color = neighbor, alpha= ICG2)) + 
  geom_point(position = position_jitter(width = 0.1, height = 0))+
  scale_color_manual(values=cbPalette[2:8], name="Neighbor")+
  facet_grid(~ anchor, switch = "x") + 
  # coord_fixed(ratio=1)+
  scale_alpha_manual(values = c(1, 0.01))+
  theme_minimal() + 
  theme(
     strip.text.x = element_text(angle = 45),   # Customize strip text
     strip.background = element_blank(),    # Style strip background
     strip.placement = "outside",                            # Place strips outside the plot area
     axis.text.x = element_blank(),
     
   )+ 
   labs(y="Directional Association Score", 
          x="Anchor")+
   guides(
     alpha = guide_legend(
       title = "Identification",
       override.aes = list(alpha =c(1, 0.2))
     ))

fig4_d
# Fig4e. boxplot of pathway ###################  

GeneList=mouseGOCC[["GOCC_GLUTAMATERGIC_SYNAPSE"]] 
dat=DSres[which(DSres$anchor=="Excitatory neuron" & DSres$neighbor=="Excitatory neuron"),]
dat$InPath=factor((dat$gene %in% GeneList), levels=c("TRUE", "FALSE"))
dat$text=ifelse(dat$DS>6, dat$gene, NA)

fig4_e=ggplot(data=dat, aes(x=InPath, y=DS)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + # Add box plot without outliers
  geom_jitter(aes( color=ICG2), position = position_jitter(width = 0.2), 
              size = 2, alpha = 0.6)  +
  scale_color_manual(values = c( "darkorange", "#999999"))+
  geom_text_repel(aes(label =text))+
  theme_classic()+
  coord_fixed(ratio=0.3)+
  labs(y="Directional Association Score", 
       x="In GO Term Glutamatergic Synapse")+
  guides(
    color = guide_legend(
      title = "Identification",
    ))
fig4_e 

# Fig4f. visualization of key gene ###################  

G=c("Cplx1", "Nsmf")
## select top 3 genes
expr=seqFISHplus_scran_sce@assays@data@listData$adjusted.counts[,seqFISHplus_scran_sce$cellClass=="Excitatory neuron"]
expr=expr[, order(colnames(expr))]
normcounts <- transform_count_to_normal(expr)

dist<-seqfish_res_list[["Excitatory neuron-Excitatory neuron"]]$distance
dist.cutoff<-seqfish_res_list[["Excitatory neuron-Excitatory neuron"]]$dist.cutoff %>% as.numeric()

## create within/outside groups
{
  gene1=data.frame(expr=as.vector(normcounts["Cplx1",]),
                   dist=factor(ifelse(dist[,2]<dist.cutoff, "Within", "Outside"),
                               levels=c("Within","Outside")))
  
  gene2=data.frame(expr=as.vector(normcounts["Nsmf",]),
                   dist=factor(ifelse(dist[,2]<dist.cutoff, "Within", "Outside"),
                               levels=c("Within","Outside")))
  
}

{
  title_size_main <- 11
  title_size <- 10
  text_size <- 10
  axis_text_size <- 9}
{
  fig4_f1 <- ggplot(gene1, aes(x = dist, y = expr)) + 
    ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, justification = -.3, point_colour = NA)+ 
    geom_boxplot( width = .25, outlier.shape = NA) + geom_point(size = 1.3, alpha = .3, position = position_jitter(seed = 1, width = .1))+
    coord_cartesian(xlim = c(1.2, NA), clip = "off") + xlab("Interaction Distance") + ylab("Expression") + 
    theme_minimal()+
    coord_fixed(ratio=0.6)+
    theme(axis.title.x = element_text(size=text_size), 
          axis.title.y = element_text(size=text_size), 
          axis.text.x = element_text(size=axis_text_size), 
          axis.text.y = element_text(size=axis_text_size), 
          plot.title = element_text(size = title_size,hjust = 0.5)) + 
    ggtitle("Cplx1")

  fig4_f2 <- ggplot(gene2, aes(x = dist, y = expr)) + 
    ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, justification = -.3, point_colour = NA)+ 
    geom_boxplot( width = .25, outlier.shape = NA) + geom_point(size = 1.3, alpha = .3, position = position_jitter(seed = 1, width = .1))+
    coord_cartesian(xlim = c(1.2, NA), clip = "off") + xlab("Interaction Distance") + ylab("Expression") + 
    theme_minimal()+
    coord_fixed(ratio=0.6)+
    theme(axis.title.x = element_text(size=text_size), 
          axis.title.y = element_text(size=text_size), 
          axis.text.x = element_text(size=axis_text_size), 
          axis.text.y = element_text(size=axis_text_size), 
          plot.title = element_text(size = title_size,hjust = 0.5)) + 
    ggtitle("Nsmf")
}
fig4_f1
fig4_f2
fig4_f= fig4_f1+ fig4_f2
fig4_f

# Output Figure3 ===============================
library(cowplot)
upper=plot_grid(fig4_a, fig4_c, labels = c("A", "C"), rel_widths = c(1,1), ncol = 2)
upper
middle=plot_grid(fig4_b, fig4_d, labels = c("B", "D"), rel_widths = c(1,3), ncol = 2)
# middle=fig4_d
middle
lower=plot_grid(fig4_e,  fig4_f, labels = c( "E",  "F"), 
                rel_widths = c(1,1.5), ncol = 2)
lower
plot_grid(upper, middle, lower,  ncol=1)


ggsave(filename = "Figures/Figure4.pdf",  # File name
       width = 11,                   # Width in inches
       height = 13,                  # Height in inches
       units = "in"                # Units: "in" (inches), "cm" (centimeters), or "mm" (millimeters)
)










