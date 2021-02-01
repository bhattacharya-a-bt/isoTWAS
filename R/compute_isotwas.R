#' Compute isoTWAS model for a set of isoforms/transcripts
#'
#' The function runs MVR models for a set of transcripts and outputs
#' the best model
#'
#' @param X matrix, design matrix of SNP dosages
#' @param Y matrix, matrix of G isoform expression across columns
#' @param Y.rep matrix, matrix of G isoform expression with replicates
#' @param R int, number of replicates
#' @param id vector, vector of sample ids showing rep to id
#' @param family character, glmnet glm family
#' @param scale logical, T/F to scale Y by Omega
#' @param alpha numeric, elastic net mixing parameter
#' @param nfolds int, number of CV folds
#' @param verbose logical
#' @param par logical, uses mclapply to parallelize model fit
#' @param n.cores int, number of parallel cores
#' @param tx_names vector, character vector of tx names in order of columns of Y
#'
#' @return optimal isoTWAS model
#'
#'
#' @export
compute_isotwas <- function(X,
                            Y,
                            Y.rep,
                            R,
                            id,
                            Omega = NULL,
                            omega_est = 'replicates',
                            omega_nlambda = 15,
                            method = c('mrce_lasso','curds_whey',
                                       'multi_enet','mrce_enet',
                                       'finemap','uni_enet',
                                       'uni_susie','uni_blup'),
                            family = 'gaussian',
                            scale = F,
                            alpha = 0.5,
                            nfolds = 5,
                            verbose,
                            par = F,
                            n.cores = NULL,
                            tx_names = NULL,
                            seed = NULL){

    ### CHECKS
    if (nrow(X) != nrow(Y)){
        stop('No. of rows of X =/= no. of rows of Y.')
    }

    N = nrow(Y)
    P = ncol(X)
    G = ncol(Y)

    if (is.null(R)){
        R = nrow(Y.full)/nrow(Y)
    } else {
        if (nrow(Y.full) != R){
            stop('No. of rows of Y.full =/= R*N')
        }
    }

    ### COMPUTE OMEGA
    omega_list = compute_omega(Y,
                               Y.rep,
                               R,
                               id,
                               method = omega_est,
                               nlambda = omega_nlambda,
                               verbose = verbose)



    return(modelList)


}
