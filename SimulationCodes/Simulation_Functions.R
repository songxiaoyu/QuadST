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

# Generate simulation data ##############################
GenSimData= function(n1, cov=F, seed=NULL){
  # fixed parameter
  k=5; R=4; rho=0.7; G=2000; nIntBlock=100;
  n=n1*k; nBlock=G/nIntBlock
  # block parameter 
  mu.b1= 0.4
  var.b1= 0.2
  a.b2= -150
  var.b2= 75
  
  # simulate spatial map
  if (cov==F) {
    x.loc=runif(n)
    y.loc=runif(n)
    celltype=paste0("CellType", sample(rep(1:k, n1), n, replace=F))
    meta=data.frame(CellID=1:n, x.loc=x.loc, y.loc=y.loc, Annotation=celltype )
    meta <-  meta %>% mutate(region = case_when( x.loc >=0 & x.loc < 0.5 & y.loc >= 0 & y.loc < 0.5 ~ "R1",
                                                 x.loc >=0.5 & x.loc <= 1 & y.loc >= 0 & y.loc < 0.5 ~ "R2",
                                                 x.loc >=0 & x.loc < 0.5 & y.loc >= 0.5 & y.loc <= 1 ~ "R3",
                                                 x.loc >=0.5 & x.loc <= 1 & y.loc >= 0.5 & y.loc <= 1 ~ "R4"))
  } else {
    nl=n+n1/2
    x.loc=runif(nl)
    y.loc=runif(nl)
    celltype=paste0("CellType", sample(rep(setdiff(1:k, 2), n1), replace=F))
    meta1=data.frame(CellID=1:nl, x.loc=x.loc, y.loc=y.loc, 
                    Annotation=c(celltype, rep("CellType2", n1*1.5)))
    meta1 <-  meta1 %>% mutate(region = case_when( x.loc >=0 & x.loc < 0.5 & y.loc >= 0 & y.loc < 0.5 ~ "R1",
                                                 x.loc >=0.5 & x.loc <= 1 & y.loc >= 0 & y.loc < 0.5 ~ "R2",
                                                 x.loc >=0 & x.loc < 0.5 & y.loc >= 0.5 & y.loc <= 1 ~ "R3",
                                                 x.loc >=0.5 & x.loc <= 1 & y.loc >= 0.5 & y.loc <= 1 ~ "R4"))
    meta2=meta1[-which(meta1$Annotation=="CellType2" & meta1$region=="R4"),]
    
    meta=rbind(meta2[meta2$Annotation!="CellType2",],
                      head(meta2[meta2$Annotation=="CellType2", ], n1))
    
  }
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
  mu1=rnorm(nIntBlock, mu.b1, var.b1) # 0.4, 0.2
  g3<- t(mvrnorm(n11, mu=mu1,Sigma=cor*2))
  g_mat[(2*nIntBlock+1):(3*nIntBlock),nn1]<-g3

  # add signal 4
  a= rnorm(nIntBlock, a.b2, var.b2) # -150, 75
  mu2=as.matrix(a) %*% (dist3[nn1]-dist.cutoff)
  g4=NULL
  for(i in 1:n11){
    tmp = mvrnorm(1, mu=mu2[,i],Sigma=cor*2)
    g4<-cbind(g4,tmp)
  }
  g_mat[(3*nIntBlock+1):(4*nIntBlock),nn1]= g4
  
  # together 
  expr[, which(meta$Annotation=="CellType1")]=g_mat
  
  if (cov==T) {
    expr[, which(meta$Annotation=="CellType1" & meta$region=="R4")]=
      expr[, which(meta$Annotation=="CellType1" & meta$region=="R4")]+0.2
  }
  return(list(meta=meta, expr=expr))
}
 
# Convert data to sce object of anchor cells 
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
RunQuadST=function(data, knn=1, bias=F,ltau=49)  {
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
  if (bias==T) {set.seed(567); sce_an$w_distance= sce_an$w_distance + 
    runif(length(sce_an$w_distance), 0, 0.005)}
  
  QRpvalue <- QuadST::test_QuadST_model(x=sce_an, 
                                        datatype="counts", 
                                        cov = NULL, 
                                        tau = dist_taus, 
                                        parallel=T)

  res <- QuadST::identify_ICGs(pMatrix=QRpvalue, fdr = 0.1)
  return(res)
}
  
## NCEM ##############################
RunNCEM <- function(data, bias=F, alpha=1, adj_method){
  meta=data$meta
  expr=data$expr
  g=nrow(expr)
  n=ncol(expr)
  y=t(expr)
  x_l= fastDummies::dummy_cols(meta$Annotation, remove_selected_columns=T) %>% as.matrix()
  # if (cov==T) {
  #   x_c= fastDummies::dummy_cols(meta$region, remove_selected_columns=T)[,-1]%>% 
  #    as.matrix() 
  # } else (x_c=NULL)
    # 
  sce_s1 <- SingleCellExperiment::SingleCellExperiment(list(counts=expr),
                                                       colData=meta,
                                                       rowData=data.frame(GeneID=rownames(expr)))
  sce_an <- QuadST::create_integrated_matrix(sce_s1, cell_id="CellID", 
                                             cell_coord1="y.loc", cell_coord2="x.loc", 
                                             cell_type="Annotation", 
                                             anchor="CellType1", neighbor= "CellType2", k=1, d.limit = Inf)
  if (bias==T) {set.seed(567); sce_an$w_distance= sce_an$w_distance + runif(length(sce_an$w_distance), 0, 0.005)}
  
  dist=sce_an$w_distance
  
  dist.cutoff=quantile(dist, probs=0.1) 
  
  dist_mat=spatstat.geom::pairdist(X=meta[, 2:3])
  diag(dist_mat)<-Inf
  A_mat=1*(dist_mat<dist.cutoff*alpha)
  x_s=1*(A_mat %*% x_l>0)
  x_ts=sapply(1:n, function(f) matrix(outer(x_l[f,], x_s[f,]))) %>% t(.)
  idx1=rep(1:5, each=5)
  idx2=rep(1:5, times=5)
  colnames(x_ts)=paste0("Type", idx1, "_to_", "Type", idx2)
  x=cbind(x_ts,x_l)
  
  ## NECM linear model regression
  res_j <- sapply(1:g, function(f) summary(lm(y[,f] ~ x-1))$coefficient["xType2_to_Type1",4])
  
  if(adj_method=="Raw"){
    return(res_j)
  }else{
    res_adj_j<-p.adjust(res_j, method=adj_method)
    return(res_adj_j)
  }
}

## Giotto ###################################
RunGiotto <- function(sce_an, alpha=1, adj_method, bias=F){
  t.test.perGene <- function(x, i, alternative = c("two.sided")){
    pvalue <- t.test(x=x[i], y = x[-i], alternative = c("two.sided"))$p.value
    return(pvalue)
  }
  w_distance=sce_an@colData@listData[["w_distance"]]
  if (bias==T) {set.seed(567); w_distance= w_distance + runif(length(w_distance), 0, 0.005)}
  
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
RunLR <- function(sce_an, adj_method, bias=F){
  y <- sce_an@colData@listData[["w_distance"]]
  
  if (bias==T) {set.seed(567); y= y + runif(length(y), 0, 0.005)}
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

# Evaluation metrics ##############################
Eva<-function(res){
  power= sum(res[1:400]<=0.1, na.rm=T)/400
  fdr = sum(res[401:4000]<=0.1, na.rm=T)/sum(res<=0.1, na.rm=T)
  # breakdown power
  power1 = sum(res[1:100]<=0.1, na.rm=T)/100
  power2 = sum(res[101:200]<=0.1, na.rm=T)/100
  power3 = sum(res[201:300]<=0.1, na.rm=T)/100
  power4 = sum(res[301:400]<=0.1, na.rm=T)/100
  return(data.frame(power=power, fdr=fdr,
                    power1=power1,power2=power2,
                    power3=power3,power4=power4))
  
}


