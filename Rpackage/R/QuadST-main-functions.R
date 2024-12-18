#' Create anchor-neighbor cell type pair integrated matrix
#`  i.e., anchor-neighbor cellpair distance and anchor cells' expression and covariates
#'
#'
#' @param x A \code{SingleCellExperiment} class.
#' @param cell_id A column name of \code{colData(object)} that stores each cell's unique id.
#' @param cell_coord1 A column name of \code{colData(object)} that stores the first spatial coordinate of each cell.
#' @param cell_coord2 A column name of \code{colData(object)} that stores the second spatial coordinate of each cell.
#' @param cell_type A column name of \code{colData(object)} that stores cell types of each cell.
#' @param anchor A name of cell type to be used for anchor cell type.
#' @param neighbor A name of cell type to be used for neighbor cell type.
#' @param k Default = 1. Find k nearest neighbors for the anchor cell.
#' @param d.limit Default=Inf. The limit of cell-cell distance for cell pairing. Cells over this distance limit
#' will not be paired even though they are the k nearest neighbors.
#' @importFrom SingleCellExperiment colData
#'
#' @return An anchor-neighbor integrated matrix in the \code{SingleCellExperiment} class.
#' @export
#'
#'
create_integrated_matrix <- function(x, cell_id, cell_coord1, cell_coord2, cell_type,
                                   anchor, neighbor, k=1, d.limit=Inf){

    object <- x
    if ( !is(object, "SingleCellExperiment") )
        stop("Object must be a SingleCellExperiment class")
    if ( !any(cell_id %in% colnames(colData(object))) )
        stop("cell_id argument must match with a column in colData(object)")
    if ( !any(cell_coord1 %in% colnames(colData(object))) )
        stop("cell_coord1 argument must match with a column in colData(object)")
    if ( !any(cell_coord1 %in% colnames(colData(object))) )
        stop("cell_coord2 argument must match with a column in colData(object)")
    if ( !any(cell_type %in% colnames(colData(object))) )
        stop("cell_type argument must match with a column in colData(object)")
    if ( !any(anchor %in% unique(colData(object)[[cell_type]])) )
        stop("anchor argument must match with an annotated cell type")
    if ( !any(neighbor %in% unique(colData(object)[[cell_type]])) )
        stop("neighbor argument must match with an annotated cell type")

    # Step 1 ---------------------------
    # Create 2d spatial point pattern object using spatstat R package
    sce_mk <- colData(object)[[cell_type]]
    sce_x <- colData(object)[[cell_coord1]]
    sce_y <- colData(object)[[cell_coord2]]
    sce_xrange <- c(min(sce_x) - 1, max(sce_x) + 1)
    sce_yrange <- c(min(sce_y) - 1, max(sce_y) + 1)
    sce_ppp <- spatstat.geom::ppp(x=sce_x, y=sce_y, xrange=sce_xrange,
                                  yrange=sce_yrange, marks=as.factor(sce_mk))
    sce_ppp[["ID"]] <- colData(object)[[cell_id]]

    # Step 2 ---------------------------
    # Find nearest source (neighbor cells) of target (anchor cells) and their distance
    nn_pairs <- .find_k_nearest_neighbors(x=sce_ppp, anchor= anchor,
                                          neighbor=neighbor, k=k, d.limit=d.limit)
    # Estimate the cell-cell interaction weighted distance
    strength=1/nn_pairs$distance
    strength_mean=tapply(strength, nn_pairs$anchor, mean)
    w_distance=1/strength_mean
    strength_k=tapply(nn_pairs$k, nn_pairs$anchor, length)
    w_distance=cbind.data.frame(anchor_idx=names(strength_mean),
                              w_distance=w_distance,
                              k=strength_k)



    # Step 3 ---------------------------
    # Subset sce object using anchor cell ids with nearest neighbor cell ids and distances.
    anchor_id <- w_distance$anchor_idx
    sce_anchor <- object[, match(anchor_id, object[[cell_id]])]
    colData(sce_anchor)[["w_distance"]] <- w_distance$w_distance
    colData(sce_anchor)[["anchor"]] <- anchor
    colData(sce_anchor)[["neighbor"]] <- neighbor

    return(sce_anchor)

}


#' Create a set of highest and lowest quantiles symmetric around median
#'
#'
#' @param min_sample_per_quantile A minimum number of samples in a given quantile level.
#' @param cell_count The number of anchor cell counts.
#' @param max_default Default = 49. The maximum number of quantile levels to use. Default = 49 means
#' quantile levels 0.02, 0.04, ..., 0.98 are under consideration.
#'
#'
#' @return A vector of quantile levels.
#' @export
#'
#'
#' @examples
#' data("seqFISHplus_scran_sce")
#' cell_id = "cellID"
#' cell_coord1 = "x"
#' cell_coord2 = "y"
#' cell_type = "cellClass"
#' anchor = "Excitatory neuron"
#' neighbor = "Astrocyte"
#' covariate = "FOV"
#' sce_an = create_cellpair_matrix(seqFISHplus_scran_sce, cell_id, cell_coord1, cell_coord2, cell_type, anchor, neighbor, cov=covariate)
#' anchor_cell_count <- length(colData(sce_an)[, cell_id])
#' dist_taus <- create_quantile_levels(min_sample_per_quantile = 5, cell_count = anchor_cell_count, max_default = 49)
#'
#'
create_quantile_levels <- function(min_sample_per_quantile, cell_count, max_default=49){

  max_ql <- max_default + 1
  number_of_quantile <- ifelse(round(cell_count/min_sample_per_quantile) < max_ql,
                               round(cell_count/min_sample_per_quantile), max_ql)
  dist_taus <- seq(0, 1, by = 1/number_of_quantile)
  dist_taus <- dist_taus[-c(1, length(dist_taus))]

  return(dist_taus)
}


#' Test anchor-neighbor distance-expression association
#'  at a set of highest and lowest quantiles symmetric around median
#'
#'
#' @param x A \code{SingleCellExperiment} class.
#' @param datatype A column name of \code{assays(object)} that stores anchor cells' gene expression levels.
#' @param cov Column names of \code{colData(object)} that needs to be adjusted as covariates.
#' @param tau A set of highest and lowest quantiles symmetric around median.
#' @param parallel Default=False. If parallel computing is activated.
#' @import QRank
#'
#' @return A matrix of quantile regression p-values in the format of genes (rows) by quantile levels (columns).
#' @export
#'
#'
test_QuadST_model <- function(x, datatype, cov=NULL, tau, parallel=F){

    object <- x
    if ( !is(object, "SingleCellExperiment") )
        stop("Object must be a SingleCellExperiment class")
    if ( !any("w_distance" %in% colnames(colData(object))) )
        stop("w_distance variable must match with a column in colData(object)")
    if ( !any(datatype %in% names(assays(object))) )
        stop("datatype argument must match with a column in colData(object)")

    # Step 1 ---------------------------
    # Set y: anchor-neighbor distance, x: anchor cells' gene expression levels, and z: covariates.
    y <- colData(object)$w_distance
    x <- t(assay(object, datatype))


    # Step 2 ---------------------------
    # Remove genes with all zeros in expression values.
    if (length(which(colSums(x) == 0)) != 0) {xMatrix <- x[, -(which(colSums(x) == 0))]
    }else{xMatrix <- x}

    if (!is.null(cov)){
        z <- colData(object)[[cov]]
        covM <- model.matrix( ~ z)[, c(-1)]
    } else {covM=NULL }

    # Step 3 ---------------------------
    # Test anchor-neighbor distance-expression association at a series of quantile levels.

    genes_wo_zeros <- colnames(xMatrix)[which(colSums(xMatrix==0)==0)]
    genes_w_zeros <- setdiff(colnames(xMatrix), genes_wo_zeros)
    pvalue=NULL
    if (length(genes_wo_zeros) != 0) {
      # Test the distance-expression association for genes with no zeros in expression values.

      if (parallel==F) {
        pvalue1 <- sapply(genes_wo_zeros, function(f)
          QRank(gene=y, snp=xMatrix[,f], cov=covM, tau=tau)$quantile.specific.pvalue) %>% t(.)
      } else {
        pvalue_temp <- parallel::mclapply(genes_wo_zeros, function(f)
          QRank(gene=y, snp=xMatrix[,f], cov=covM, tau=tau)$quantile.specific.pvalue)
        pvalue1=unlist(pvalue_temp) %>% matrix(., ncol=length(tau), byrow = T)
        rownames(pvalue1)=genes_wo_zeros
      }
      pvalue=rbind(pvalue, pvalue1)

    }
    if (length(genes_w_zeros)!=0) {
        # Test distance-expression association for genes with some zeros in expression values.
      if (parallel==F) {
        pvalue2 <- sapply(genes_w_zeros, function(f) .QRank_multi(y=y, x=cbind(xMatrix[,f], 1*I(xMatrix[,f] != 0)),
                                                                  cov=covM, tau=tau, alternative="two-sided-directional")$quantile.specific.pvalue) %>% t(.)

      } else {
        pvalue_temp <- parallel::mclapply(genes_w_zeros, function(f) .QRank_multi(y=y, x=cbind(xMatrix[,f], 1*I(xMatrix[,f] != 0)),
                                                                  cov=covM, tau=tau, alternative="two-sided-directional")$quantile.specific.pvalue)
        pvalue2=unlist(pvalue_temp) %>% matrix(., ncol=length(tau), byrow = T)
        rownames(pvalue2)=genes_w_zeros
      }
      pvalue=rbind(pvalue, pvalue2)

    }
    # Combine p-values for genes with no zeros and with some zeros in expression values.
    colnames(pvalue)=tau
    return(pvalue)
}


#' Identify cell-cell interaction changes genes (ICGs)
#'
#'
#' @param pMatrix A matrix of p-values in the genes (rows) by quantile levels (columns) format.
#' @param fdr Default= 0.1. The nominal false discovery rate to claim significance.
#'
#'
#' @return a list of inference results
#' \item{summary.table}{A summary of the results, including the quantile index and level of the
#' cell-cell interaction, and the No. of ICGs significant at this quantile level.}
#' \item{data.table}{A table of gene-level results, including gene name, p-value at the specified quantile
#' level, the empirical FDR, and whether it's a ICG or not.}
#' @export
#'
#'
identify_ICGs <- function(pMatrix, fdr = 0.1){

  pvalue <- pMatrix
  if ( !is(pvalue, "matrix") ) stop("Object must be a matrix")

  L <- ncol(pvalue)
  M <- floor(L/2)
  eFDR =NoSig = NULL
  for (m in 1:M) {
    cuts=pvalue[,m]

    pB=pvalue[,c(1:m)] # near (large signal- ICG)
    pA=pvalue[,seq(L, L-m+1)] # further (little signal)

    nA=sapply(cuts, function(f) mean(pA<f))
    nB=sapply(cuts, function(f) mean(pB<f))
    eFDR1= (nA+0.000001)/(nB+0.000001)

    o <- order(cuts, decreasing = F)
    ro <- order(o)
    eFDR1_clean=pmin(1, cummin(eFDR1[o]))[ro]
    names(eFDR1_clean)=names(eFDR1)
    eFDR=cbind(eFDR, eFDR1_clean)
    NoSig=c(NoSig, sum(eFDR1_clean<fdr))


    # eFDR1_reorder=eFDR1[order(eFDR1)]
    # o <- order(cuts)
    # ro <- order(o)
    # eFDR1_clean=pmin(1, eFDR1_reorder)[ro]
    # names(eFDR1_clean)=names(eFDR1)
    # temp=data.frame(cuts, nA, nB, eFDR1, eFDR1_clean)

    eFDR=cbind(eFDR, eFDR1_clean)
    NoSig=c(NoSig, sum(eFDR1_clean<fdr))
  }

  if (any(NoSig>0)) {
    m1=which.max(NoSig)
    s.table=data.frame(idx_ICG=m1,
                       Q_taus=colnames(pvalue)[m1],
                       sig_gene_count=NoSig[m1])
    d.table=data.frame(gene=rownames(pvalue),
                       pvalue=pvalue[,s.table$Q_taus],
                       eFDR=eFDR[,m1],
                       ICG=1*(eFDR[,m1]<fdr))

  }else{
    m1=which.min(apply(eFDR, 2, min))
    s.table =  data.frame(idx_ICG=NA,
                          Q_taus=colnames(pvalue)[m1],
                          sig_gene_count=0)

    d.table =  data.frame(gene=rownames(pvalue),
                          pvalue=pvalue[,s.table$Q_taus],
                          eFDR=eFDR[,m1],
                          ICG=0)
  }


  return(list(summary.table=s.table, data.table=d.table))

}

#' Identify cell-cell interaction distance
#'
#'
#' @param x A \code{SingleCellExperiment} class.
#' @param ICG.summary ICG summary.
#'@param k K nearest neighbor.
#'
#' @return Interacting distance.
#' @export
#'
#'
ICG_distance <- function(x, ICG.summary, k) {
  object <- x
  if ( !is(object, "SingleCellExperiment") ) stop("Object must be a SingleCellExperiment class")
  if (k>1)  stop("Distance can only be identified if 1 nearest neighor is used")

  # Step 1 ---------------------------
  # Set y: anchor-neigbhor distance, x: anchor cells'gene expression levels, and z: covariates.
  y <- object$w_distance
  dist=quantile(y, probs=as.numeric(ICG.summary$Q_taus))
  return(dist)

}
