rm(list=ls())
library(tidyverse)
library(ggplot2)
library(patchwork)
library(gridExtra)
setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/SpatialTranscriptomics/Paper_QuadST_Revision")

set.colors <- c("#31A354", "#3182BD", "#F0E442", "#d7191c")

## Fig2a #########################

dat=read.csv("Results/Fig2a_Simulation_Summary.csv")
qt=c("0.25, 0.5, 0.75",
     "0.2, 0.4, ..., 0.8",
     "0.1, 0.2, ..., 0.9", 
     "0.05, 0.1, ..., 0.95", 
     "0.02, 0.04, ..., 0.98", 
     "0.01, 0.02, ..., 0.99", 
     "0.001, 0.002, ..., 0.999")
dat=dat %>%   mutate(Power=round(mean_power*100,1), 
                     FDR= round(mean_fdr*100,1)) %>%
  mutate(Quantile=factor(qt, levels=qt))

pa1=ggplot(data=dat, aes(x=Quantile, y=FDR, fill=Quantile))+
  geom_bar(stat = "identity")+
  theme_classic() + 
  lims(y=c(0,100))+ 
  geom_hline(yintercept=10, linetype="dashed")+
  geom_text(aes(label = FDR), vjust = -0.5)+
  facet_wrap(~"FDR") +
  theme(plot.title = element_text(hjust = 0.5, size=12),
        legend.position="none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y=element_blank(),
        strip.background = element_rect(fill = "gray80"))


pa2=ggplot(data=dat, aes(x=Quantile, y=Power, fill=Quantile))+
  geom_bar(stat = "identity")+
  geom_text(aes(label = Power), vjust = -0.5)+
  theme_classic() + 
  lims(y=c(0,100))+ 
  facet_wrap(~"Power") +
  theme(legend.position="none",
        axis.title.y=element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "gray80"))
ppa1=pa1+ggtitle("Quantile Level")
fig_a = ppa1/pa2
fig_a
## Fig2b:  #########################

dat=read.csv("Results/Fig2b_Simulation_Summary.csv") %>%
  mutate(Power=round(mean_power*100,1), 
         FDR= round(mean_fdr*100,1))

pb1=ggplot(data=dat, aes(x=knn, y=FDR, fill=knn))+
  geom_bar(stat = "identity")+
  theme_classic() + 
  lims(y=c(0,100))+ 
  geom_hline(yintercept=10, linetype="dashed")+
  geom_text(aes(label = FDR), vjust = -0.5)+
  facet_wrap(~"FDR") +
  theme(plot.title = element_text(hjust = 0.5, size=12),
        legend.position="none",
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y=element_blank(),
        strip.background = element_rect(fill = "gray80"))

pb2=ggplot(data=dat, aes(x=knn, y=Power, fill=knn))+
  geom_bar(stat = "identity")+
  geom_text(aes(label = Power), vjust = -0.5)+
  theme_classic() + 
  lims(y=c(0,100))+ 
  facet_wrap(~"Power") +
  theme(legend.position="none",
        axis.title.y=element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = "gray80"))

ppb1=pb1+ggtitle("KNN Selection")
fig_b=ppb1/pb2
fig_b

## Fig2c #########################

dat=read.csv("Results/Fig2c_Simulation_Summary.csv") %>%
  mutate(Power=round(mean_power*100,1), 
         FDR= round(mean_fdr*100,1),
         n=as.factor(n))
pc1=ggplot(data=dat, aes(x=n, y=FDR, fill=n))+
  geom_bar(stat = "identity")+
  theme_classic() + 
  lims(y=c(0,100))+ 
  geom_hline(yintercept=10, linetype="dashed")+
  geom_text(aes(label = FDR), vjust = -0.5)+
  facet_wrap(~"FDR") +
  scale_fill_brewer(palette="Dark2")+ 
  theme(plot.title = element_text(hjust = 0.5, size=12),
        legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y=element_blank(),
        strip.background = element_rect(fill = "gray80"))

pc2=ggplot(data=dat, aes(x=n, y=Power, fill=n))+
  geom_bar(stat = "identity")+
  geom_text(aes(label = Power), vjust = -0.5)+
  theme_classic() + 
  lims(y=c(0,100))+ 
  facet_wrap(~"Power") +
  scale_fill_brewer(palette="Dark2")+ 
  theme(legend.position="none",
        axis.title.y=element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = "gray80"))

ppc1=pc1+ggtitle("Sample Size")
fig_c=ppc1/pc2
fig_c
## Fig2d #########################
dat=read.csv("Results/Fig2d_Simulation_Summary.csv") 

dat=dat%>% 
  mutate(Method = factor(method, levels=c("Giotto", "NCEM", "LR", "QuadST")), 
         Power=round(mean_power*100,1), 
         FDR= round(mean_fdr*100,1))

p1=ggplot(data=dat, aes(x=Method, y=FDR, fill=Method))+
  geom_bar(stat = "identity")+
  theme_classic() + 
  lims(y=c(0,100))+ 
  scale_fill_manual(values=set.colors) + 
  geom_hline(yintercept=10, linetype="dashed")+
  geom_text(aes(label = FDR), vjust = -0.5)+
  facet_wrap(~"FDR") +
  theme(plot.title = element_text(hjust = 0.5, size=12),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.y=element_blank(),
        strip.background = element_rect(fill = "gray80"))

p2=ggplot(data=dat, aes(x=Method, y=Power, fill=Method))+
  geom_bar(stat = "identity")+
  geom_text(aes(label = Power), vjust = -0.5)+
  theme_classic() + 
  lims(y=c(0,100))+ 
  scale_fill_manual(values=set.colors)+ 
  facet_wrap(~"Power") +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.background = element_rect(fill = "gray80"))


pp1=p1+ggtitle("Ideal Case")
fig_d=pp1/p2

## Fig2e #########################
dat=read.csv("Results/Fig2e_Simulation_Summary.csv")  

dat=dat%>% 
  mutate(Method = factor(method, levels=c("Giotto", "NCEM", "LR", "QuadST")), 
         Power=round(mean_power*100,1), 
         FDR= round(mean_fdr*100,1))

p1=ggplot(data=dat, aes(x=Method, y=FDR, fill=Method))+
  geom_bar(stat = "identity")+
  theme_classic() + 
  lims(y=c(0,100))+ 
  scale_fill_manual(values=set.colors) + 
  geom_hline(yintercept=10, linetype="dashed")+
  geom_text(aes(label = FDR), vjust = -0.5)+
  facet_wrap(~"FDR") +
  theme(plot.title = element_text(hjust = 0.5, size=12),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.y=element_blank(),
        strip.background = element_rect(fill = "gray80"))

p2=ggplot(data=dat, aes(x=Method, y=Power, fill=Method))+
  geom_bar(stat = "identity")+
  geom_text(aes(label = Power), vjust = -0.5)+
  theme_classic() + 
  lims(y=c(0,100))+ 
  scale_fill_manual(values=set.colors)+ 
  facet_wrap(~"Power") +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.background = element_rect(fill = "gray80"))


pp1=p1+ggtitle("Inaccurate Neigbhor Estimate")
fig_e=pp1/p2

## Fig2f: Measurement Error in Distance #########################
dat=read.csv("Results/Fig2f_Simulation_Summary.csv")  

dat=dat%>% 
  mutate(Method = factor(method, levels=c("Giotto", "NCEM", "LR", "QuadST")), 
         Power=round(mean_power*100,1), 
         FDR= round(mean_fdr*100,1))

p1=ggplot(data=dat, aes(x=Method, y=FDR, fill=Method))+
  geom_bar(stat = "identity")+
  theme_classic() + 
  lims(y=c(0,100))+ 
  scale_fill_manual(values=set.colors) + 
  geom_hline(yintercept=10, linetype="dashed")+
  geom_text(aes(label = FDR), vjust = -0.5)+
  facet_wrap(~"FDR") +
  theme(plot.title = element_text(hjust = 0.5, size=12),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.y=element_blank(),
        strip.background = element_rect(fill = "gray80"))

p2=ggplot(data=dat, aes(x=Method, y=Power, fill=Method))+
  geom_bar(stat = "identity")+
  geom_text(aes(label = Power), vjust = -0.5)+
  theme_classic() + 
  lims(y=c(0,100))+ 
  scale_fill_manual(values=set.colors)+ 
  facet_wrap(~"Power") +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.background = element_rect(fill = "gray80"))

pp1=p1+ggtitle("Measurement Error in Distance")
fig_f=pp1/p2
## Fig2g: Under Unmeasured Confounder #########################
dat=read.csv("Results/Fig2g_Simulation_Summary.csv") 

dat=dat%>% 
  mutate(Method = factor(method, levels=c("Giotto", "NCEM", "LR", "QuadST")), 
         Power=round(mean_power*100,1), 
         FDR= round(mean_fdr*100,1))

p1=ggplot(data=dat, aes(x=Method, y=FDR, fill=Method))+
  geom_bar(stat = "identity")+
  theme_classic() + 
  lims(y=c(0,100))+ 
  scale_fill_manual(values=set.colors) + 
  geom_hline(yintercept=10, linetype="dashed")+
  geom_text(aes(label = FDR), vjust = -0.5)+
  facet_wrap(~"FDR") +
  theme(plot.title = element_text(hjust = 0.5, size=12),
        legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        axis.title.y=element_blank(),
        strip.background = element_rect(fill = "gray80"))

p2=ggplot(data=dat, aes(x=Method, y=Power, fill=Method))+
  geom_bar(stat = "identity")+
  geom_text(aes(label = Power), vjust = -0.5)+
  theme_classic() + 
  lims(y=c(0,100))+ 
  scale_fill_manual(values=set.colors)+ 
    facet_wrap(~"Power") +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        strip.background = element_rect(fill = "gray80"))

pp1=p1+ggtitle("Under Unmeasured Confounder")
fig_g=pp1/p2
## Putting together #########################

library(cowplot)
upper=plot_grid(fig_a, fig_b, fig_c, labels = c("A", "B", "C"), ncol = 3)
lower=plot_grid(fig_d, fig_e, fig_f, fig_g,  labels = c( "D", "E", "F", "G"), ncol = 4)
plot_grid(upper, lower,  ncol=1)

ggsave(filename = "Figures/Figure2.pdf",  # File name
  width = 11,                   # Width in inches
  height = 8,                  # Height in inches
  units = "in"                # Units: "in" (inches), "cm" (centimeters), or "mm" (millimeters)
)
