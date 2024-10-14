#' Find nearest neighbor cells and their distances from anchor cells
#'
#'
#' @param x A point pattern object of class 'ppp'.
#' @param source_cell A name of cell type for neighbor cells.
#' @param target_cell A name of cell type for anchor cells.
#' @importFrom magrittr "%>%"
#' @return A dataframe with anchor-neighbor cell pair with distance.
#'
#'
#' @noRd
#'
#'
.find_nearest_neighbors <- function(x, source_cell, target_cell){

    object <- x
    if ( !is(object, "ppp") )
        stop("Object must be a ppp class: representing a point pattern dataset in the two-dimensional plane.")

    # Find the nearest neighbor cell (source_cell) and its distance for each anchor cell (target_cell)
    nn <- spatstat.geom::nnwhich(object, by=as.factor(object$marks), k=1)
    nnd <- spatstat.geom::nndist(object, by=as.factor(object$marks), k=1)
    tcell <- which(object$marks == target_cell)
    nn_source <- nn[tcell, source_cell]
    nnd_source <- nnd[tcell, source_cell]
    nn_pair <- cbind.data.frame(target=tcell, source=nn_source, distance=nnd_source) %>% dplyr::arrange(distance)

    return(nn_pair)
}


#' Calculate test statistics of the parameters in QuadST quantile regression model: Signed-QRank
#'
#'
#' @param y Anchor-neighbor cell-pair distance.
#' @param x Anchor cells' gene expression values.
#' @param cov Covariates to adjust for.
#' @param tau A set of quantile levels at which test statistics are calculated.
#' @param alternative By default, use "two-sided-directional".
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


#' Create a set of highest and lowest quantiles symmetric around median
#'
#'
#' @param min_sample_per_quantile A minimum number of samples in a given quantile level.
#' @param cell_count The number of anchor cell counts.
#' @param max_default The maximum number of quantile levels to use by dafault.
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
create_quantile_levels <- function(min_sample_per_quantile, cell_count, max_default){

    max_ql <- max_default + 1
    number_of_quantile <- ifelse(round(cell_count/min_sample_per_quantile) < max_ql,
                                 round(cell_count/min_sample_per_quantile), max_ql)
    dist_taus <- seq(0, 1, by = 1/number_of_quantile)
    dist_taus <- dist_taus[-c(1, length(dist_taus))]

    return(dist_taus)
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


#' Combine gene-specific p-values at each highest and lowest quantiles symmetrically around the median quantile (0.5)
#'  using Cauchy combination test (ACAT R package)
#'
#'
#' @param x A quantile regression p-value matrix: genes by quantiles.
#' @import ACAT
#'
#' @return A Cauchy combination test p-value matrix: genes by quantiles (without the median quantile).
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
#' https://stackoverflow.com/questions/16357962/r-split-numeric-vector-at-position/16358095#16358095
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


#' Identify ICGs with an empirically estimated FDR p-value threshold
#'
#'
#' @param x A Caucy combination test p-value matrix: genes by quantiles (without the median quantile)
#' @param pvalue_cutoff An initial p-value threshold value. Use 0.05 by default. This serves as a starting point to find a significant p-value with the empirical FDR procedure.
#' @param ABratio_cutoff A norminal false discovery rate to control. Use 0.1 by default.
#' @param ABconst A constant used in the empirical FDR procedure. Use 0.1 by default.
#'
#'
#' @return a list of results
#' \item{ABratio}{}
#' \item{p_cutoff}{}
#' \item{sig_gene_count}{}
#' \item{q_status}{}
#' \item{sig_gene}{}
#' \item{sig_gene_id}{}
#'
#'
#' @noRd
#'
#'
.control_eFDR <- function(x, fdr = fdr){

    pvalue <- x
    if ( !is(pvalue, "matrix") ) stop("Object must be a matrix")

    L <- ncol(pvalue)
    M <- floor(L/2)
    eFDR =NoSig = vector("list", length = M)
    for (m in 1:M) {
      pB=pvalue[,m] # near - ICG
      pA=pvalue[,(L-m)] # further
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
      eFDR[[m]]=eFDR1_clean
      NoSig[[m]]=sum(eFDR1_clean<ABratio_cutoff)
    }

    if (any(unlist(NoSig)>0)) {
      idx_ICG=which.max(NoSig)
      Q_taus=colnames(pvalue)[idx_ICG]
      pB_ICG=pvalue[,idx_ICG]
      eFDR_ICG=eFDR[[idx_ICG]]
      sig_gene_count=sum(eFDR_ICG<fdr)
      sig_gene_id=names(which(eFDR_ICG<fdr))
      sig_gene_data=data.frame(Gene=sig_gene_id, pvalue=pB_ICG[which(eFDR_ICG<fdr)],
                               eFDR=eFDR_ICG[which(eFDR_ICG<fdr)])
      sig_gene_data=sig_gene_data[order(sig_gene_data$pvalue),]

    }else{

      Q_taus =      sig_gene_data =      sig_gene_count =      sig_gene_id = NA
    }

    q_index <- 1:M

    return(list(q_status=Q_taus, sig_gene_count=sig_gene_count,
                sig_gene_data=sig_gene_data, sig_gene_id=sig_gene_id))
}


#' Calculate empirically estimated FDR (A to B ratio) between each highest and lowest quantiles symmetrically around the median.
#'
#'
#'
#' @param pvalue A Caucy combination test p-value matrix: genes by quantiles.
#' @param pvalue A Caucy combination test p-value matrix: genes by quantiles.
#' @param ABconst A constant used in empirically estimated FDR calculation.
#'
#'
#' @return Empirically estimated FDR (A to B ratio) across quantiles
#'
#'
#' @noRd
#'
#'
.cal_ABratio <- function(x, p_cutoff, ABconst = 0.1){

    pvalue <- x
    if ( !is(pvalue, "matrix") )
     stop("Object must be a matrix")

    # Step 1 ---------------------------
    # Transform Caucy combination test p-value matrix: genes by quantiles (without the median quantile)
    # into a matrix with its element indicating p-value is lower than a pvalue cutoff (p_cutoff).
    L <- ncol(pvalue)
    q <- round(L/2)
    pvalue_sig <- I(pvalue <= p_cutoff) * 1

    # Step 2 ---------------------------
    # Calculate empirically estimated FDR (A to B ratio where A and B refer to higher and lower quantiles respectively).
    if (q == 1){
        countA <- sum(pvalue_sig[, L:(q+1)])
        countB <- sum(pvalue_sig[, 1:q])
        ABratio <- (countA + ABconst)/(countB + ABconst)
        ABratio <- is.finite(ABratio)*ABratio
        names(ABratio) <- colnames(pvalue)[1:q]
    }else if (q > 1){
        countA <- colSums(pvalue_sig[, L:(q+1)])
        countB <- colSums(pvalue_sig[, 1:q])
        ABratio <- (countA + ABconst)/(countB + ABconst)
        ABratio <- is.finite(ABratio)*ABratio
        names(ABratio) <-  colnames(pvalue)[1:q]
    }

    return(ABratio)
}
