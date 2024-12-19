rm(list=ls())
# suppressMessages(library(QuadST))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(parallel))
suppressMessages(library(dplyr))
suppressMessages(library(autoReg))
suppressMessages(library(MASS))
suppressMessages(library(spatstat))
suppressMessages(library(rlist))
suppressMessages(library(mvtnorm))
setwd("/Users/songxiaoyu152/NUS Dropbox/Xiaoyu Song/SpatialTranscriptomics/Paper_QuadST_Revision")
source("../Paper_QuadST/Github/Rpackage/R/QuadST-helper-functions.R")
source("../Paper_QuadST/Github/Rpackage/R/QuadST-main-functions.R")
source("../Github/SimulationCodes/Simulation_Functions.R")

## Method comparison between QuadST and others (Figure 2d-g) ##############################
#  Fig 2d: ideal case
# {n1=1000; run_knn=1; run_ltau=49; alpha=1; cov=F; bias=F; mc=10} 
# # Fig 2e: inaccurate neigbhor estimate (note: one can run NCEM and Giotto only)
# {n1=1000;run_knn=1; run_ltau=49; cov=F; alpha=1.4; bias=F; mc=10} 
# # Fig 2f: biased cell-cell distance
{n1=1000;run_knn=1; run_ltau=49; cov=F; alpha=1; bias=T; mc=10} 
# # Fig 2g: unmeasured covariates
# {n1=1000;run_knn=1; run_ltau=49; cov=T; alpha=1; bias=F; mc=10} 


set.seed(123)
seed=round(runif(mc)*10000)
dat_list=NULL
for (i in 1:mc) {
  dat_i=GenSimData(n1 = n1, seed=seed[i], cov=cov)
  dat_list[[i]]=dat_i
}
sce_an_list=Gen_sce_an_list(dat_list)



QuadST_out=NULL
for (i in 1:mc) {
  cat(i,"\n")
  dat_i=dat_list[[i]]
  res_i=RunQuadST(data=dat_i, knn=run_knn, ltau=run_ltau, bias=bias)  
  out_i=data.frame(method="QuadST",n=n1, knn=run_knn, 
                   ltau=run_ltau, Eva(res_i$data.table$eFDR)) 
  QuadST_out=rbind(QuadST_out, out_i)
}
QuadST_out

NCEM_out=NULL
for (i in 1:mc) {
  dat_i=dat_list[[i]]
  cat(i,"\n")
  res_i=RunNCEM(data=dat_i, alpha=alpha, adj_method="BH", bias=bias)  
  out_i=data.frame(method="NCEM",alpha=alpha, adj_method="BH",Eva(res_i)) 
  NCEM_out=rbind(NCEM_out, out_i)
}
NCEM_out
QuadST_out

Giotto_out=NULL
for (i in 1:mc) {
  sce_an_i=sce_an_list[[i]]
  res_i=RunGiotto(sce_an=sce_an_i, alpha=alpha, adj_method="BH", bias=bias)  
  out_i=data.frame(method="Giotto",alpha=alpha, adj_method="BH",Eva(res_i)) 
  Giotto_out=rbind(Giotto_out, out_i)
}
Giotto_out

LR_out=NULL
for (i in 1:mc) {
  sce_an_i=sce_an_list[[i]]
  res_i=RunLR(sce_an=sce_an_i, adj_method="BH", bias=bias)  
  out_i=data.frame(method="LR",adj_method="BH",Eva(res_i)) 
  LR_out=rbind(LR_out, out_i)
}
LR_out

sel_cols=c("method","power","fdr") 
raw_out=rbind(QuadST_out[,sel_cols],NCEM_out[,sel_cols], 
              Giotto_out[,sel_cols], LR_out[,sel_cols])
#outDir=file.path("Results/Fig2d_Simulation_Results.csv")
# outDir=file.path("Results/Fig2e_Simulation_Results.csv")
# outDir=file.path("Results/Fig2f_Simulation_Results.csv")
outDir=file.path("Results/Fig2g_Simulation_Results.csv")
write.table(raw_out, file=outDir, col.names = TRUE, row.names = FALSE, sep = ",")
# raw_out=read.csv("Results/Fig2d_Simulation_Results.csv")
library(dplyr)
summary_out = raw_out %>% group_by(method) %>% 
  summarize(mean_power = mean(power), 
            mean_fdr = mean(fdr))
# outDir=file.path("Results/Fig2d_Simulation_Summary.csv")
# outDir=file.path("Results/Fig2e_Simulation_Summary.csv")
# outDir=file.path("Results/Fig2f_Simulation_Summary.csv")
outDir=file.path("Results/Fig2g_Simulation_Summary.csv")
write.table(summary_out, file=outDir, col.names = TRUE, row.names = FALSE, sep = ",")
summary_out


## Evalation of QuadST  (Figure 2a-c)  ##############################

#  Fig 2a: No. of taus
{n1=1000; run_knn=1; run_ltau=c(3, 4, 9, 19, 49, 99, 999); mc=10}
set.seed(234)
seed=round(runif(mc)*10000)
dat_list=NULL
for (i in 1:mc) {
  dat_i=GenSimData(n1 = 1000, seed=seed[i], cov=cov)
  dat_list[[i]]=dat_i
}
QuadST_ltau=NULL
for (l in run_ltau) {
  for (i in 1:mc) {
    dat_i=dat_list[[i]]
    res_i=RunQuadST(data=dat_i, knn=run_knn, ltau=l)  
    out_i=data.frame(method="QuadST",n=n1, mc=i, knn=run_knn, 
                     ltau=l, Eva(res_i$data.table$eFDR)) 
    QuadST_ltau=rbind(QuadST_ltau, out_i)
  }
}

outDir=file.path("Results/Fig2a_Simulation_Results.csv")
write.table(QuadST_ltau, file=outDir, col.names = TRUE, row.names = FALSE, sep = ",")

summary_out = QuadST_ltau %>% group_by(ltau) %>% 
  summarize(mean_power = mean(power), 
            mean_fdr = mean(fdr))
outDir=file.path("Results/Fig2a_Simulation_Summary.csv")
write.table(summary_out, file=outDir, col.names = TRUE, row.names = FALSE, sep = ",")
summary_out

#  Fig 2b: No. of taus
{n1=1000; run_knn=c(1,2,3,4,5); run_ltau=49; mc=10}

set.seed(234)
seed=round(runif(mc)*10000)
dat_list=NULL
for (i in 1:mc) {
  dat_i=GenSimData(n1 = 1000, seed=seed[i])
  dat_list[[i]]=dat_i
}
sce_an_list=Gen_sce_an_list(dat_list)

QuadST_knn=NULL
for (k in run_knn) {
  print(c("k:", k))
  for (i in 1:mc){
    print(c("i:", i))
    dat_i=dat_list[[i]]
    res_i=RunQuadST(data=dat_i, knn=k, ltau=run_ltau)  
    out_i=data.frame(method="QuadST", mc=i, n=n1, knn=k, 
                     ltau=run_ltau, Eva(res_i$data.table$eFDR)) 
    QuadST_knn=rbind(QuadST_knn, out_i)
  }
}

outDir=file.path("Results/Fig2b_Simulation_Results.csv")
write.table(QuadST_knn, file=outDir, col.names = TRUE, row.names = FALSE, sep = ",")

summary_out = QuadST_knn %>% group_by(knn) %>% 
  summarize(mean_power = mean(power), 
            mean_fdr = mean(fdr))
outDir=file.path("Results/Fig2b_Simulation_Summary.csv")
write.table(summary_out, file=outDir, col.names = TRUE, row.names = FALSE, sep = ",")
summary_out

#  Fig 2c: sample size
{n1=c(100, 300, 500, 700, 1000); run_knn=1; run_ltau=49; mc=10}

QuadST_n=NULL
for (nn in n1) {
  print(c("nn:", nn))
  set.seed(234)
  seed=round(runif(mc)*10000)
  for (i in 1:mc) {
    print(c("i:", i))
    dat_i=GenSimData(n1 = nn, seed=seed)
    res_i=RunQuadST(data=dat_i, knn=run_knn, ltau=run_ltau)  
    out_i=data.frame(method="QuadST",n=nn, mc=i, knn=run_knn, 
                     ltau=run_ltau, Eva(res_i$data.table$eFDR)) 
    QuadST_n=rbind(QuadST_n, out_i)
  }
}

outDir=file.path("Results/Fig2c_Simulation_Results.csv")
write.table(QuadST_n, file=outDir, col.names = TRUE, row.names = FALSE, sep = ",")

summary_out = QuadST_n %>% group_by(n) %>% 
  summarize(mean_power = mean(power), 
            mean_fdr = mean(fdr))
outDir=file.path("Results/Fig2c_Simulation_Summary.csv")
write.table(summary_out, file=outDir, col.names = TRUE, row.names = FALSE, sep = ",")
summary_out










