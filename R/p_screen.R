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
p_screen <- function(Z,
                     Sigma){

    level1=3e4
    cutoff2=.05
    level2=3e5
    cutoff3=.01
    level3=3e6
    cutoff4=1e-5
    level4=3e7

    t = length(Z)
    pvals = 2*pnorm(-abs(Z))
    minp=min(pvals)
    rank=order(pvals)
    ordered=pvals[rank]
    v = cov2cor(Sigma)

    lower=rep(qnorm(minp/2),t)
    upper=rep(qnorm(1-minp/2),t)

    p_ACT= 1 - mvtnorm::pmvnorm(lower=lower,
                                upper=upper,
                                sigma=v,
                                maxpts=level1,
                                abseps=1e-13)
    if (p_ACT < cutoff2) {
        p_ACT = 1 - mvtnorm::pmvnorm(lower=lower,
                                     upper=upper,
                                     sigma=v,
                                     maxpts=level2,
                                     abseps=1e-13)
        if (p_ACT < cutoff3) {
            p_ACT = 1 - mvtnorm::pmvnorm(lower=lower,
                                         upper=upper,
                                         sigma=v,
                                         maxpts=level3,
                                         abseps=1e-13)
            if (p_ACT < cutoff4) {
                p_ACT = 1 - mvtnorm::pmvnorm(lower=lower,
                                             upper=upper,
                                             sigma=v,
                                             maxpts=level4,
                                             abseps=1e-13)
            }
        }
    }

    return(p_ACT[1])

}
