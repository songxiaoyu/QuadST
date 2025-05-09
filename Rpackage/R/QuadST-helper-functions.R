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
#' @export
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
      nn_pair1=cbind.data.frame(anchor=object$ID[acell], neighbor=nn1,
                                distance=nnd1, k=i)
      nn_pair=rbind(nn_pair, nn_pair1)
    }
    nn_pair2=nn_pair[which(nn_pair$distance < d.limit),]

    return(nn_pair2)
}


#' Signed-QRank to calculate test statistics of multiple parameters in QuadST quantile regression model
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
#' @export
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




#' Transform expression data to a normally distribution with the minimum value at zero
#'
#'
#' @param x An G by N expression matrix.
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
    N=ncol(x)
    G=nrow(x)
    # Transform a given distribution to a normal distribution.
    x_norm <- qnorm(1:N/(N+1))
    expr=matrix(NA, ncol=N, nrow=G)
    for (i in 1:G) {
      q1=rank(x[i,])
      x2=x_norm[q1]
      min=min(x2)
      expr[i,]=x2-min
    }
    colnames(expr)=colnames(x)
    rownames(expr)=rownames(x)
    return(expr)
}


