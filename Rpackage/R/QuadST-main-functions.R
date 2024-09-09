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
#' @param cov Column names of \code{colData(object)} that need to be adjusted as covariates.
#'
#'
#' @return A \code{SingleCellExperiment} class with an anchor-neighbor cell type pair integrated matrix
#' @export
#'
#'
create_cellpair_matrix <- function(x, cell_id, cell_coord1, cell_coord2, cell_type,
                                   anchor, neighbor, cov=NULL){

    object <- x
    if ( !is(object, "SingleCellExperiment") )
        stop("Object must be a SingleCellExperiment class")
    if ( !any(cell_id %in% colnames(SingleCellExperiment::colData(object))) )
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
    if ( !all(cov %in% colnames(colData(object))) )
        stop("cov argument must match with columns in colData(object)")

    # Step 1 ---------------------------
    # Create 2d spatial point pattern object using spatstat R package
    sce_mk <- colData(object)[[cell_type]]
    sce_x <- colData(object)[[cell_coord1]]
    sce_y <- colData(object)[[cell_coord2]]
    sce_xrange <- c(min(sce_x) - 1, max(sce_x) + 1)
    sce_yrange <- c(min(sce_y) - 1, max(sce_y) + 1)
    sce_ppp <- spatstat.geom::ppp(x=sce_x, y=sce_y, xrange=sce_xrange, yrange=sce_yrange, marks=sce_mk)
    sce_ppp[[cell_id]] <- colData(object)[[cell_id]]
    sce_ppp[[cov]] <- colData(object)[[cov]]
    # References for incorporating 3d spatial point pattern
    # https://inside.mines.edu/~jdzimmer/tutorials/Section1.html
    # https://github.com/aproudian2/rapt/tree/master/R

    # Step 2 ---------------------------
    # Find nearest source (neighbor cells) of target (anchor cells) and their distance
    source_cell <- neighbor
    target_cell <- anchor
    nn_pairs <- .find_nearest_neighbors(sce_ppp, source_cell, target_cell)
    distance <- nn_pairs$distance
    anchor_id <- colData(object)[[cell_id]][nn_pairs$target]

    # Step 3 ---------------------------
    # Subset sce object using anchor cell ids with nearest neigbhor cell ids and distances.
    sce_anchor <- object[, object[[cell_id]] %in% anchor_id]
    sce_anchor <- sce_anchor[, match(anchor_id, sce_anchor[[cell_id]])]
    colData(sce_anchor)$distance <- distance
    colData(sce_anchor)$anchor <- anchor
    colData(sce_anchor)$neighbor <- neighbor

    return(sce_anchor)
    # Use the following if spatial point patten of cells needs to be returend.
    #return(list(sce_an=sce_anchor, sce_ppp=sce_ppp))
}


#' Test anchor-neighbor distance-expression association
#'  at a set of highest and lowest quantiles symmetric around median
#'
#'
#' @param x A \code{SingleCellExperiment} class.
#' @param dist A column name of \code{colData(object)} that stores anchor-neighbor cell pair distances.
#' @param expr A column name of \code{assays(object)} that stores anchor cells' gene expression levels.
#' @param cov Column names of \code{colData(object)} that needs to be adjusted as covariates.
#' @param tau A set of highest and lowest quantiles symmetric around median.
#' @import QRank
#'
#' @return A matrix of quntile regression p-values: genes (rows) by quantiles (columns).
#' @export
#'
#'
test_QuadST_model <- function(x, dist, expr, cov=NULL, tau){

    object <- x
    if ( !is(object, "SingleCellExperiment") )
        stop("Object must be a SingleCellExperiment class")
    if ( !any(dist %in% colnames(colData(object))) )
        stop("dist argument must match with a column in colData(object)")
    if ( !any(expr %in% names(assays(object))) )
        stop("expr argument must match with a column in colData(object)")

    # Step 1 ---------------------------
    # Set y: anchor-neigbhor distance, x: anchor cells'gene expression levels, and z: covariates.
    y <- colData(object)[[dist]]
    x <- t(assay(object, expr))
    z <- colData(object)[[cov]]

    # Step 2 ---------------------------
    # Remove genes with all zeros in expression values.
    if (length(which(colSums(x) == 0)) != 0) {
        xMatrix <- x[, -(which(colSums(x) == 0))]
    }else{
        xMatrix <- x
    }
    if (!is.null(cov)){
        covM <- model.matrix( ~ z)[, c(-1)]
    }

    # Step 3 ---------------------------
    # Test anchor-neighbor distance-expression association at a series of quantile levels.
    if (length(which(colSums(xMatrix==0)==0)) != 0) {
        # Test the distance-expression association for genes with no zeros in expression values.
        genes_wo_zeros <- colnames(xMatrix)[which(colSums(xMatrix==0)==0)]
        pvalue <- tryCatch(sapply(genes_wo_zeros, function(f) QRank(gene=y, snp=xMatrix[,f], cov=covM, tau=tau)$quantile.specific.pvalue) %>% t(.), error = function(e) NULL)

        # Test the distance-expression association for genes with some zeros in expression values.
        genes_w_zeros <- setdiff(colnames(xMatrix), genes_wo_zeros)
        pvalue1 <- sapply(genes_w_zeros, function(f) tryCatch(.QRank_multi(y=y, x=cbind(xMatrix[,f], 1*I(xMatrix[,f] != 0)), cov=covM, tau=tau, alternative="two-sided-directional")$quantile.specific.pvalue, error = function(e) NULL), simplify=FALSE)
        genes_w_zeros <- names(pvalue1)[!sapply(pvalue1, is.null)]
        pvalue1 <- as.matrix(dplyr::bind_cols(pvalue1[genes_w_zeros])) %>% t()

        # Combine p-values for genes with no zeros and with some zeros in expression values.
        genes_w_QRpvalue <- colnames(xMatrix)[colnames(xMatrix) %in% c(genes_wo_zeros, genes_w_zeros)]
        pvalue <- rbind(pvalue, pvalue1)[genes_w_QRpvalue,]
    }else{
        # Test distance-expression association for genes with some zeros in expression values.
        genes_w_zeros <- colnames(xMatrix)
        pvalue <- sapply(genes_w_zeros, function(f)
          tryCatch(.QRank_multi(y=y, x=cbind(xMatrix[,f], 1*I(xMatrix[,f] != 0)),
                                cov=covM, tau=tau, alternative="two-sided-directional")$quantile.specific.pvalue,
                   error = function(e) NULL), simplify=FALSE)
        genes_w_zeros <- names(pvalue)[!sapply(pvalue, is.null)]
        pvalue <- as.matrix(dplyr::bind_cols(pvalue[genes_w_zeros])) %>% t()
        pvalue <- matrix(pvalue, ncol=length(tau), dimnames=list(genes_w_zeros, tau))
    }

    return(pvalue)
}


#' Identify cell-cell interaction changes genes (IGGs)
#'
#'
#' @param x A \code{SingleCellExperiment} class.
#' @param y A matrix of quntile regression p-values: genes (rows) by quantiles (columns).
#' @param dist A column name of \code{colData(object1)} that stores anchor-neighbor cell pair distances.
#' @param expr A column name of \code{assays(object1)} that stores anchor cells' gene expression levels.
#' @param cov Column names of \code{colData(object1)} that need to be adjusted as covariates.
#' @param tau A set of quantile levels at which test statistics are calculated.
#' @param p_thres An initial p-value threshold value. Use 0.05 by default. This serves as a starting point to find a significant p-value with the empirical FDR procedure.
#' @param fdr A norminal false discovery rate to control. Use 0.1 by default.
#' @param ABconst A constant used in the empirical FDR procedure. Use 0.1 by default.
#'
#'
#' @return a list of inference results
#' \item{ACAT}{Combined gene-specific p-values at each lowest and highest quantile symmetrically around the median quantile using Cauchy Combination test.}
#' \item{p_sig}{Significant p value threshold across quantile levels at a prescribed FDR threshold, e.g., 0.1.}
#' \item{q_sig}{Logicals that indicate if there exist significant p value thresholds across quantile levels.}
#' \item{ICGs}{Identified cell-cell interaction changed gene names.}
#' \item{q_int}{An interaction quantile.}
#' \item{dist_int}{An interaction distance.}
#' \item{DA_score}{Directional association scores calculated at an interaction quantile.}
#' @export
#'
#'
identify_ICGs <- function(x, y, dist, expr, cov = NULL, tau, p_thres = 0.05, fdr = 0.1, ABconst = 0.1){

    object1 <- x
    object2 <- y
    if ( !is(object1, "SingleCellExperiment") )
        stop("Object must be a SingleCellExperiment class")
    if ( !any(dist %in% colnames(colData(object1))) )
        stop("dist argument must match with a column in colData(object)")
    if ( !any(expr %in% names(assays(object1))) )
        stop("expr argument must match with a column in colData(object)")
    if ( !all(cov %in% colnames(colData(object1))) )
        stop("cov argument must match with columns in colData(object)")
    if ( !is(object2, "matrix") )
     stop("Object must be a matrix")

    # Step 1 ---------------------------
    # Calculate combined gene-specific p-values across each highest and lowest quantiles symmetrically around
    # the median quantile (0.5).
    ACATpvalue <- .ACAT_QRpvalue(object2)

    # Step 2 ---------------------------
    # Identify ICGs controlling empirical false discovery rate
    res <- .control_eFDR(ACATpvalue, pvalue_cutoff = p_thres, ABratio_cutoff = fdr, ABconst = ABconst)
    ABratio <- res[["ABratio"]]
    q_status <- res[["q_status"]]
    q_sig <- sapply(q_status, function(x) !is.na(x))
    p_cutoff <- res[["p_cutoff"]]
    sig_gene_count <- res[["sig_gene_count"]]
    sig_gene_pvalue <- res[["sig_gene"]]
    sig_gene_id <- res[["sig_gene_id"]]
    q_int <- which.max(sig_gene_count)

    # Step 3 ---------------------------
    # Calculate directional association scores if an cell-cell interaction quantile (q_int) is identified.
    if (any(q_sig)){
        # Calculate an anchor-neighbor cell interaction distance from the identified cell-cell interaction quantile (q_int).
        dist_int <- quantile(colData(object1)[, dist], prob=as.numeric(names(q_int)))

        # Calculate directional association scores at the identified cell-cell interaction quantile (q_int).
        y <- colData(object1)[[dist]]
        x <- t(assay(object1, expr))
        z <- colData(object1)[[cov]]
        covM <- model.matrix( ~ z)[,-c(1)]
        q_tau <- as.numeric(names(q_int))

        genes <- rownames(object1)
        tstat <- sapply(genes, function(f) tryCatch(.QRank_multi(y=y, x=cbind(x[,f], 1*I(x[,f] != 0)), cov=covM, tau=q_tau, alternative="two-sided-directional")$quantile.specific.test["x",], error = function(e) NULL), simplify=TRUE)
        q_ACATpvalue <- ACATpvalue[genes,  as.character(q_tau)]
        q_tstat <- unlist(tstat)[names(q_ACATpvalue)]
        DAscore <- -log10(q_ACATpvalue)*sign(q_tstat)
    }else{
        dist_int <- NA
        DAscore <- NA
    }

    return(list(ACAT=ACATpvalue, p_sig=p_cutoff, q_sig=q_sig, ICGs=sig_gene_id, q_int=q_int, dist_int=dist_int, DA_score=DAscore))
}
