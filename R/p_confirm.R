#' Compute screening P-value for one gene
#'
#' The function runs a modified version of P_ACT to test
#' hypothesis of at least one gene-transcript association
#'
#' @param Z numeric, vector of TWAS Z-scores
#' @param Sigma numeric, matrix of TWAS LD within a gene
#'
#' @return screen P-value for one gene
#'
#'
#' @importFrom mvtnorm pmvnorm
#'
#'
#' @export
p_confirm <- function(p,
                      alpha){

    o <- order(p)
    n <- length(p)
    if (n==1) {
        adjustment=0
    } else {
        adjustment=c(n-1,(n-1):1)
    }
    pAdjusted <- p[o]*adjustment
    pAdjusted <- pmin(pAdjusted,1)
    pAdjusted <- cummax(pAdjusted)
    pBack <- vector(length=length(p))
    pBack[o] <- pAdjusted
    names(pBack) <- names(p)
    pBack[pBack > alpha] = 1
    return(pBack)

}


