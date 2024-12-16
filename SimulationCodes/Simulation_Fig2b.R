rm(list=ls())
# suppressMessages(library(QuadST))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(parallel))
suppressMessages(library(dplyr))
suppressMessages(library(autoReg))
suppressMessages(library(MASS))
suppressMessages(library(spatstat))
suppressMessages(library(rlist))
mainDir <- "/Users/yuqingshang/Dropbox/Paper_QuadST_Revision"
source("../Paper_QuadST/Github/Rpackage/R/QuadST-helper-functions.R")
source("../Paper_QuadST/Github/Rpackage/R/QuadST-main-functions.R")

# Generate simulation data ##############################
GenSimData= function(n1 = 1000, seed=NULL){
  # fixed parameter
  k=5; R=4; rho=0.7; G=4000; nBlock=40
  n=n1*k; nIntBlock=G/nBlock
  # block parameter 
  mu.b1= 0.4
  var.b1= 0.2
  a.b2= -150
  var.b2= 75
  
  # simulate spatial map
  x.loc=runif(n)
  y.loc=runif(n)
  celltype=paste0("CellType", sample(rep(1:k, n1), n, replace=F))
  
  meta=data.frame(CellID=1:n, x.loc=x.loc, y.loc=y.loc, Annotation=celltype )
  meta <-  meta %>% mutate(region = case_when( x.loc >=0 & x.loc < 0.5 & y.loc >= 0 & y.loc < 0.5 ~ "R1",
                                               x.loc >=0.5 & x.loc <= 1 & y.loc >= 0 & y.loc < 0.5 ~ "R2",
                                               x.loc >=0 & x.loc < 0.5 & y.loc >= 0.5 & y.loc <= 1 ~ "R3",
                                               x.loc >=0.5 & x.loc <= 1 & y.loc >= 0.5 & y.loc <= 1 ~ "R4"))
 
  # simulate expression profile 
  expr=matrix(NA, G, n)
  colnames(expr)=paste0("Cell", meta$CellID)
  rownames(expr)=paste0("Gene", 1:G)
  
  # null data for all cell types 
  autoRegCor = function(sqmatrix, rho) {# auto correlation parameter
    dim=dim(sqmatrix)[1]
    for (i in 1:dim) {
      for (j in 1:dim) {
        sqmatrix[i,j]=sqmatrix[j,i]=rho^(abs(i-j))
      }
    }
    return(sqmatrix)
  }
  cor= autoRegCor(sqmatrix=diag(nIntBlock), rho=rho)
  for (j in 1:k) { # k cell types
    g0=NULL
    for(i in 1:nBlock){
      t=t(mvrnorm(n1, mu=rep(0, nIntBlock), cor))
      g0=rbind(g0, t)
    }
    expr[, which(meta$Annotation==paste0("CellType", j))]=g0
  }
    
  #  add signal ------- 
  # Cal nearest cell pairs.
  dist=spatstat.geom::pairdist(X=meta[, 2:3]) # dist: cal cell distance for all cell types
  dist2=dist[which(meta$Annotation=="CellType1"), 
             which(meta$Annotation=="CellType2")] # dist2: select cell dist of anchor and neighbor cell types
  # determine CCI mechanisms 1 -->  1-nearest neighbor within 10% of cells.  
  dist3=apply(dist2, 1, min) # dist3: nearest neighbor distance for each anchor cell
  dist.cutoff=quantile(dist3, probs=0.1)
  nn1=which(dist3<dist.cutoff) # nn1: idx in meta1 for anchor cells having interactions
  n11=length(nn1) # n11: number of interactions
  
  dist4=apply(dist2, 1, function(f) sort(f)[2])
  dist.cutoff2=quantile(dist4, probs=0.1)
  nn2=which(dist4<dist.cutoff2) # nn2: idx in meta1 for anchor cells having interactions
  n12=length(nn2) # n12: number of interactions
  
  # update 
  g_mat<-expr[, which(meta$Annotation==paste0("CellType", 1))]

  # add signal 1 -- Having Cell 2 in the neighbor impact the expression of cell 1 
  mu1=rnorm(nIntBlock, mu.b1, var.b1) # 10, 0.2
  g1<- t(mvrnorm(n11, mu=mu1,Sigma=cor)) 
  g_mat[1:nIntBlock, nn1]<-g1
  
  # add signal 2
  a= rnorm(nIntBlock, a.b2, var.b2) # -150, 75
  mu2=as.matrix(a) %*% (dist3[nn1]-dist.cutoff)
  g2=NULL
  for(i in 1:n11){
    tmp = mvrnorm(1, mu=mu2[,i],Sigma=cor)
    g2<-cbind(g2,tmp)
  }
  g_mat[(nIntBlock+1):(2*nIntBlock),nn1]<- g2
    
  # add signal 3
  mu1=rnorm(nIntBlock, mu.b1, var.b1) # 10, 0.2
  g3<- t(mvrnorm(n11, mu=mu1,Sigma=cor*4))
  g_mat[(2*nIntBlock+1):(3*nIntBlock),nn1]<-g3
  
  # add signal 4
  a= rnorm(nIntBlock, a.b2, var.b2) # -150, 75
  mu2=as.matrix(a) %*% (dist3[nn1]-dist.cutoff)
  g4=NULL
  for(i in 1:n11){
    tmp = mvrnorm(1, mu=mu2[,i],Sigma=cor*4)
    g4<-cbind(g4,tmp)
  }
  g_mat[(3*nIntBlock+1):(4*nIntBlock),nn1]= g4
  
  # add signal 5
  a= rnorm(nIntBlock, a.b2, var.b2) # -150, 75
  mu2=as.matrix(a) %*% (dist4[nn2]-dist.cutoff2)
  g5=NULL
  for(i in 1:n12){
    tmp = mvrnorm(1, mu=mu2[,i],Sigma=cor)
    g5<-cbind(g5,tmp)
  }
  g_mat[(4*nIntBlock+1):(5*nIntBlock),nn2]= g5
  
  
  # together 
  expr[, which(meta$Annotation=="CellType1")]=g_mat
  
  return(list(meta=meta, expr=expr))
}
 
# Convert data to sce object of anchor cells ##############################

Gen_sce_an_list<-function(dat_list){
  sce_an_list=NULL
  for(i in 1: length(dat_list)){
    dat_i<-dat_list[[i]]
    cell_feature=dat_i$meta
    expr=dat_i$expr
    sce_s1 <- SingleCellExperiment::SingleCellExperiment(list(counts=expr),
                                                         colData=cell_feature,
                                                         rowData=data.frame(GeneID=rownames(expr)),
                                                         metadata="Setting1")
    sce_an <- QuadST::create_integrated_matrix(sce_s1, cell_id="CellID", 
                                               cell_coord1="y.loc", cell_coord2="x.loc", 
                                               cell_type="Annotation", 
                                               anchor="CellType1", neighbor= "CellType2", k=1, d.limit = Inf)
    sce_an_list[[i]]=sce_an
  }
  return(sce_an_list)
}
# QuadST ##############################
RunQuadST=function(data, knn=1, ltau=49)  {
  # Load Data------------------
  cell_feature=data$meta
  expr=data$expr
 
  dist_taus=1:ltau/(1+ltau)
  # create SingleCellExperiment object 
  sce_s1 <- SingleCellExperiment::SingleCellExperiment(list(counts=expr),
                                 colData=cell_feature,
                                 rowData=data.frame(GeneID=rownames(expr)),
                                 metadata="Setting1")
  sce_an <- QuadST::create_integrated_matrix(sce_s1, cell_id="CellID", 
                                             cell_coord1="y.loc", cell_coord2="x.loc", 
                                             cell_type="Annotation", 
                                     anchor="CellType1", neighbor= "CellType2", k=knn, 
                                     d.limit = Inf)
  QRpvalue <- QuadST::test_QuadST_model(x=sce_an, datatype="counts", 
                                cov = NULL, tau = dist_taus, parallel=T)
  
  res <- QuadST::identify_ICGs(pMatrix=QRpvalue, fdr = 0.1)
  return(res)
}
  
## NCEM ##############################
RunNCEM <- function(data, alpha=1, adj_method){
  meta=data$meta
  expr=data$expr
  g=nrow(expr)
  n=ncol(expr)
  y=t(expr)
  x_l= fastDummies::dummy_cols(meta$Annotation, remove_selected_columns=T) %>% as.matrix()
  # 
  sce_s1 <- SingleCellExperiment::SingleCellExperiment(list(counts=expr),
                                                       colData=meta,
                                                       rowData=data.frame(GeneID=rownames(expr)))
  sce_an <- QuadST::create_integrated_matrix(sce_s1, cell_id="CellID", 
                                             cell_coord1="y.loc", cell_coord2="x.loc", 
                                             cell_type="Annotation", 
                                             anchor="CellType1", neighbor= "CellType2", k=1, d.limit = Inf)
  dist=sce_an$w_distance
  dist.cutoff=quantile(dist, probs=0.1) 
  dist_mat=spatstat.geom::pairdist(X=meta[, 2:3])
  A_mat=1*(dist_mat<dist.cutoff*alpha)
  x_s=1*(A_mat %*% x_l>0)
  x_ts=sapply(1:n, function(f) matrix(outer(x_l[f,], x_s[f,]))) %>% t(.)
  idx1=rep(1:5, each=5)
  idx2=rep(1:5, times=5)
  colnames(x_ts)=paste0("Type", idx1, "_to_", "Type", idx2)
  x=cbind(x_l, x_ts)
  ## NECM linear model regression
  res_j <- sapply(1:g, function(f) summary(lm(y[,f] ~ x-1))$coefficient[10,4])
  
  if(adj_method=="Raw"){
    return(res_j)
  }else{
    res_adj_j<-p.adjust(res_j, method=adj_method)
    return(res_adj_j)
  }
}

## Giotto ###################################
RunGiotto <- function(sce_an, alpha=1, adj_method){
  t.test.perGene <- function(x, i, alternative = c("two.sided")){
    pvalue <- t.test(x=x[i], y = x[-i], alternative = c("two.sided"))$p.value
    return(pvalue)
  }
  w_distance=sce_an@colData@listData[["w_distance"]]
  dist.cutoff=quantile(w_distance, probs=0.1)
  
  y <- sce_an@colData@listData[["w_distance"]]
  expr <- sce_an@assays@data@listData[["counts"]]
  Interact_cells <- which(y < alpha*dist.cutoff)
  ## Giotto T-test
  res_j <- apply(expr, 1, function(gene) 
    t.test.perGene(gene, i=Interact_cells, alternative =c("two.sided")) )
  if(adj_method=="Raw"){
    return(res_j)
  }else{
    res_adj_j<-p.adjust(res_j, method=adj_method)
    return(res_adj_j)
  }
}

## LR ###################################
RunLR <- function(sce_an, adj_method){
  y <- sce_an@colData@listData[["w_distance"]]
  expr <- sce_an@assays@data@listData[["counts"]]
  k=nrow(expr)
  ## Linear Regression
  res_j <- sapply(1:k, function(f) summary(lm(y ~ expr[f,]))$coefficient[2,4])
  if(adj_method=="Raw"){
    return(res_j)
  }else{
    res_adj_j<-p.adjust(res_j, method=adj_method)
    return(res_adj_j)
  }
}

# # Comparison ##############################
# Eva<-function(res){
#   power= sum(res[1:400]<=0.1)/400
#   fdr = sum(res[401:4000]<=0.1, na.rm=T)/sum(res<=0.1, na.rm=T)
#   # breakdown power
#   power1 = sum(res[1:100]<=0.1)/100
#   power2 = sum(res[101:200]<=0.1)/100
#   power3 = sum(res[201:300]<=0.1)/100
#   power4 = sum(res[301:400]<=0.1)/100
#   return(data.frame(power=power, fdr=fdr,
#                     power1=power1,power2=power2,
#                     power3=power3,power4=power4))
#   #return(data.frame(power=power, fdr=fdr))
# }

# Comparison ##############################
Eva<-function(res){
  power= sum(res[1:500]<=0.1)/500
  fdr = sum(res[501:4000]<=0.1, na.rm=T)/sum(res<=0.1, na.rm=T)
  # breakdown power
  power1 = sum(res[1:100]<=0.1)/100
  power2 = sum(res[101:200]<=0.1)/100
  power3 = sum(res[201:300]<=0.1)/100
  power4 = sum(res[301:400]<=0.1)/100
  power5 = sum(res[401:500]<=0.1)/100
  return(data.frame(power=power, fdr=fdr,
                    power1=power1,power2=power2,
                    power3=power3,power4=power4,
                    power5=power5))
  #return(data.frame(power=power, fdr=fdr))
}


## Use function (1) data generation  ##############################

## params #############################
{run_knn=1; run_ltau=49;mc=1}


set.seed(123)
seed=round(runif(mc)*10000)
dat_list=NULL
for (i in 1:mc) {
  dat_i=GenSimData(n1 = 1000, seed=seed[i])
  dat_list[[i]]=dat_i
}
sce_an_list=Gen_sce_an_list(dat_list)

## run QuadST ##############################
quadst_list=QuadST_out=NULL
for (i in 1:mc) {
  cat(i,"\n")
  dat_i=dat_list[[i]]
  res_i=RunQuadST(data=dat_i, knn=1, ltau=99)  
  quadst_list[[i]]=res_i
  out_i=data.frame(method="QuadST",n=1000, knn=run_knn, 
                   ltau=run_ltau, Eva(res_i$data.table$eFDR)) 
  QuadST_out=rbind(QuadST_out, out_i)
}


NCEM_out=NULL
for (i in 1:mc) {
  sce_an_i=sce_an_list[[i]]
  cat(i,"\n")
  res_i=RunNCEM(data=dat_i, alpha=1, adj_method="BH")  
  out_i=data.frame(method="NCEM_A",alpha=1, adj_method="BH",Eva(res_i)) 
  NCEM_out=rbind(NCEM_out, out_i)
  
  #res_i=RunNCEM(data=dat_i, alpha=1.4, adj_method="BH")
  #out_i=data.frame(method="NCEM_B",alpha=1.4, adj_method="BH",Eva(res_i))
  #NCEM_out=rbind(NCEM_out, out_i)
}
QuadST_out
NCEM_out


outDir=file.path(mainDir,"Data/simulation/quadst_10_var2.RData")
save(quadst_list,file=outDir)
outDir=file.path(mainDir,"Figures/tmp_quadst.csv")
write.table(QuadST_out, file=outDir,col.names = TRUE, row.names = FALSE, sep = ",")

tmp_by_block=NCEM_out %>%group_by(method) %>% summarise(across(power:power4, mean))

tmp_by_block=QuadST_out %>%group_by(method) %>% summarise(across(power:power4, mean))


Giotto_out=NULL
for (i in 1:mc) {
  sce_an_i=sce_an_list[[i]]
  
  res_i=RunGiotto(sce_an=sce_an_i, alpha=1, adj_method="BH")  
  out_i=data.frame(method="Giotto_A",alpha=1, adj_method="BH",Eva(res_i)) 
  Giotto_out=rbind(Giotto_out, out_i)
  
  res_i=RunGiotto(sce_an=sce_an_i, alpha=1.4, adj_method="BH")  
  out_i=data.frame(method="Giotto_B",alpha=1.4, adj_method="BH",Eva(res_i)) 
  Giotto_out=rbind(Giotto_out, out_i)
}

LR_out=NULL
for (i in 1:mc) {
  sce_an_i=sce_an_list[[i]]
  res_i=RunLR(sce_an=sce_an_i, adj_method="BH")  
  out_i=data.frame(method="LR",adj_method="BH",Eva(res_i)) 
  LR_out=rbind(LR_out, out_i)
}

## summarize output #########################
sel_cols=c("method","power","fdr") 
raw_out=rbind(QuadST_out[,sel_cols],NCEM_out[,sel_cols], Giotto_out[,sel_cols], LR_out[,sel_cols])
outDir=file.path(mainDir,"Figures/raw_compare.csv")
write.table(raw_out, file=outDir,col.names = TRUE, row.names = FALSE, sep = ",")

library(dplyr)
summary_out = raw_out %>% group_by(method) %>% summarize(mean_power = mean(power), mean_fdr = mean(fdr))
outDir=file.path(mainDir,"Figures/summary_compare.csv")
write.table(summary_out, file=outDir,col.names = TRUE, row.names = FALSE, sep = ",")

