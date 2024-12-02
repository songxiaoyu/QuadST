#' Find k nearest neighbor cells and their distances from anchor cells
#'
#'
#' @param x A point pattern object of class 'ppp'.
#' @param anchor A name of cell type for anchor cells.
#' @param neighbor A name of cell type for neighbor cells.
#' @param k Default = 1. Find k nearest neighbors for the anchor cell.
#' @param d.limit Default=Inf. The limit of cell-cell distance for cell pairing. Cells over this distance limit
#' will not be paired even though they are the k nearest neighbors.
#' @importFrom magrittr "%>%"
#' @return A data frame for anchor-neighbor cell pairs with their distances and interaction strengths.
#'
#'
#' @noRd
#'
#'
.find_k_nearest_neighbors <- function(x, anchor, neighbor, k, d.limit){

    object <- x
    if ( !is(object, "ppp") )
        stop("Object must be a ppp class: representing a point pattern dataset in the two-dimensional plane.")

    # Find the k nearest neighbor cell within a distance limit for each anchor cell
    facet=as.factor(object$marks==neighbor)
    acell <- which(object$marks == anchor)
    nn_pair=NULL
    for (i in 1:k) {
      nn1 <- spatstat.geom::nnwhich(object, by=facet, k=i)[acell,2]
      nnd1 <- spatstat.geom::nndist(object, by=facet, k=i)[acell,2]
      nn_pair1=cbind.data.frame(anchor=object[[cell_id]][acell], neighbor=nn1,
                                distance=nnd1, k=i)
      nn_pair=rbind(nn_pair, nn_pair1)
    }
    nn_pair2=nn_pair[which(nn_pair$distance < d.limit),]

    return(nn_pair2)
}


#' Calculate test statistics of the parameters in QuadST quantile regression model: Signed-QRank
#'
#'
#' @param y An variable indicating anchor-neighbor interaction strength.
#' @param x  Expression level of a gene in anchor cell type.
#' @param cov Covariates to adjust for.
#' @param tau A set of quantile levels at which test statistics are calculated.
#' @param alternative Default ="two-sided-directional". Possible values include "two-sided-no-direction",
#' "two-sided-directional", "greater-directional", and "less-directional", providing the flexibility on
#' the treatment of directions and one-sided vs two-sided tests.
#'
#' @import QRank quantreg
#' @return A list of Signed-QRank test results
#' \item{quantile.specific.pvalue}{Signed-QRank p-value at each quantile.}
#' \item{quantile.specific.test}{Signed-QRank test statistics at each quantile.}
#'
#'
#' @noRd
#'
#'
.QRank_multi=function (y, x, cov = NULL, tau,
                      alternative=c("two-sided-no-direction",
                                    "two-sided-directional",
                                    "greater-directional",
                                    "less-directional")) {

  ltau = length(tau)
  x = as.matrix(x)
  y = as.matrix(y)
  zz = cbind(rep(1, nrow(y)), cov)
  p=ncol(x)

  xstar = lm(x ~ zz - 1)$residual
  Q=t(xstar) %*% xstar

  # Calculate p-values and test statistics
  ltau=length(tau)
  pvalue=NULL
  tstat=NULL
  for (l in 1:ltau) {
    tau1=tau[l]
    ranks = suppressWarnings(quantreg::rq.fit.br(zz, y, tau = tau1)$dual -  (1 - tau1))
    Sn = as.matrix(t(xstar) %*% (ranks))
    VN = tau1*(1-tau1) * Q

    if (alternative=="two-sided-no-direction") {
      test=t(Sn) %*% solve(VN) %*% Sn
      pvalue=c(pvalue, pchisq(test, df=p, lower.tail = F))
    } else {
      z_each=Sn/sqrt(diag(VN))
      r=length(Sn)
      rho=cov2cor(VN)[lower.tri(VN)]
      z_bar=sum(Sn/sqrt(diag(VN)))/sqrt(r+2*sum(rho))
      if (alternative=="two-sided-directional") {
        pvalue=c(pvalue, 2*pnorm(abs(z_bar), lower.tail = F))
      }
      if (alternative=="greater") {
        pvalue=c(pvalue, pnorm(z_bar, lower.tail = F))
      }
      if (alternative=="less") {
        pvalue=c(pvalue, pnorm(z_bar, lower.tail = T))
      }
      tstat=cbind(tstat, c(z_each,  z_bar))
    }
  }
  names(pvalue)=tau
  colnames(tstat)=tau
  rownames(tstat)=c(paste0("x", 1:length(z_each)), "x")

  return(list(quantile.specific.pvalue = pvalue, quantile.specific.test = tstat))
}




#' Transform cell-specific bias adjusted counts to a normally distribution with the minimum value at zero
#'
#'
#' @param x An expression matrix.
#'
#'
#' @return A transformed expression matrix.
#' @export
#'
#'
#' @seealso [test_QuadST_model()] to see how this function is used.
#'
#'
transform_count_to_normal <- function(x){

    # Step 1 ---------------------------
    # Transform a given distribution to a normal distribution.
    fun_ecdf <- ecdf(x)
    x_ecdf <- fun_ecdf(x)
    x_norm <- qnorm(x_ecdf)

    # Step 2 ---------------------------
    # If the maximum of ecdf value (x_ecdf) is equal to 1, the inverse of cdf (x_norm) can become infinity.
    # Ensure that the maximum of ecdf value (x_ecdf) is less than 1.
    if ( any(!is.finite(x_norm)) ){
        x_num <- c(as.numeric(x), max(as.numeric(x) + 1))
        fun_ecdf <- ecdf(x_num)
        x_ecdf <- fun_ecdf(x_num)[-length(as.numeric(x_num))]
        x_norm <- qnorm(x_ecdf)
    }

    # Step 3 ---------------------------
    # Set the minimum value of a normal distribution at 0.
    x_norm_minzero <- matrix(x_norm, nrow=nrow(x), ncol=ncol(x), dimnames=list(rownames(x), colnames(x))) - min(x_norm)

    return(x_norm_minzero)
}


#' Smoothing p-values
#'
#' Combine gene-specific p-values at each highest and lowest quantiles symmetrically around the median quantile (0.5)
#' using Cauchy combination test (ACAT R package)
#'
#'
#' @param x A quantile regression p-value matrix: genes by quantiles.
#' @import ACAT
#'
#' @return A Cauchy combination of the p-value matrix: genes by quantile levels (without the median quantile).
#'
#'
#' @noRd
#'
#'
.ACAT_QRpvalue <- function(x){

        object <- x
        if ( !is(object, "matrix") )
         stop("Object must be a matrix.")
        if ( rlang::is_empty(colnames(object)) )
         stop("Matrix column names should specificy the quantile level of each column.")

        y_taus <- colnames(object)
        L <- length(y_taus)

        # If the length of quantile levels is even:
        if (L %% 2 == 1) {
            qC <- ceiling(L/2)
            qLR <- .splitAt(1:L, qC)
            qL <- qLR[[1]]
            qR <- setdiff(qLR[[2]], qC)
        }else{
            qC <- L/2 + 1
            qLR <- .splitAt(1:L, qC)
            qL <- qLR[[1]]
            qR <- qLR[[2]]
        }

        # If the length of quantile levels is odd:
        if (length(qL) == 1) {
            QR_LR <- object[,c(qL, qR)]
        }else if (length(qL) > 1) {
            QR_L <- apply(object, 1, function(q) sapply(qL, function(z) ACAT::ACAT(q[1:z])))
            rownames(QR_L) <- y_taus[qL]
            QR_L <- QR_L %>% t(.)
            QR_R <- apply(object, 1, function(q) sapply(qR, function(z) ACAT::ACAT(q[z:L])))
            rownames(QR_R) <- y_taus[qR]
            QR_R <- QR_R %>% t(.)
            QR_LR <- cbind(QR_L, QR_R)
        }

        return(QR_LR)
}


#' Splitting a vector into two at a certain index
#'
#'
#' @param x A numeric vector
#' @param pos An index to split a vector into two sub vectors.
#'
#'
#' @return A list consisting of two vectors
#'
#'
#' @noRd
#'
#'
.splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))


#' Identify ICGs using empirical FDR.
#'
#'
#' @param x A smoothed p-value matrix in the format of genes by quantile levels (without the median quantile)
#' @param fdr Default =0.1. A nominal false discovery rate to control.
#'
#'
#' @return a list of results
#' \item{summary.table}{A summary of the results, including the quantile index and level of the
#' cell-cell interaction, and the No. of ICGs significant at this quantile level.}
#' \item{data.table}{A table of gene-level results, including gene name, p-value at the specified quantile
#' level, the empirical FDR, and whether it's a ICG or not.}
#'
#'
#' @noRd
#'
#'
.control_eFDR <- function(x, fdr=0.1 ){

    pvalue <- x
    if ( !is(pvalue, "matrix") ) stop("Object must be a matrix")

    L <- ncol(pvalue)
    M <- floor(L/2)
    eFDR =NoSig = NULL
    for (m in 1:M) {
      pB=pvalue[,m] # near (large signal- ICG)
      pA=pvalue[,(L-m+1) ] # further (little signal)

      cuts=pB*1.0001
      # cuts=cuts[which(cuts<pvalue_cutoff)]
      # cuts=cuts[order(cuts)]
      nA=sapply(cuts, function(f) sum(pA<f))
      nB=sapply(cuts, function(f) sum(pB<f))
      eFDR1= (nA+0.001)/(nB+0.001)

      o <- order(pB, decreasing = TRUE)
      ro <- order(o)
      eFDR1_clean=pmin(1, cummin(eFDR1[o]))[ro]
      names(eFDR1_clean)=names(eFDR1)
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


