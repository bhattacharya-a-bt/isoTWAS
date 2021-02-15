#' Multivariate SuSiE model
#'
#' The function trains multivariate SuSiE for all
#' isoform transcripts jointly
#'
#' @param X matrix, design matrix of SNP dosages
#' @param Y matrix, matrix of G isoform expression across columns
#' @param nfolds int, number of CV folds
#' @param verbose logical
#' @param par logical, uses mclapply to parallelize model fit
#' @param n.cores int, number of parallel cores
#' @param tx_name vector, character vector of tx names in order of columns of Y
#' @param seed int, random seed
#'
#' @return data frame of elastic net, lasso, and LMM based predictions
#'
#' @importFrom mmbr msusie
#' @importFrom pbapply pbapply
#' @importFrom tibble tibble
#' @importFrom rlist list.append
#'
#' @export
multivariate_mmbr <- function(X,
                              Y,
                              nfolds = 5,
                              verbose = T,
                              tx_names = NULL,
                              L = 10,
                              seed){

    if (!is.null(colnames(Y))){
        tx_names = colnames(Y)
    }


    Y = as.data.frame(Y)
    if (!is.null(tx_names)){
        colnames(Y) = tx_names
    }

    set.seed(seed)
    train.folds = caret::createFolds(1:nrow(Y),
                                     k = nfolds,
                                     returnTrain = T)
    pred = matrix(ncol = ncol(Y),
                  nrow = nrow(Y))
    for (tr in 1:nfolds){
        Y.tr = Y[train.folds[[tr]],]
        X.tr = X[train.folds[[tr]],]
        X.test = X[-train.folds[[tr]],]
        sink(tempfile())
        prior_covar = mmbr::create_mash_prior(sample_data =
                                                  list(X=X.tr,
                                                       Y=Y.tr,
                                                       residual_variance =
                                                           cov(Y.tr)),
                                              max_mixture_len=-1)
        m <- mmbr::msusie(X = X.tr,
                          Y = Y.tr,
                          prior_variance = prior_covar,
                          precompute_covariances = T,
                          L = 10,
                          intercept = T,
                          standardize = T)
        sink()
        pred[-train.folds[[tr]],] <- mmbr::predict.mmbr(m, X.test)
    }
    r2.vec = sapply(1:ncol(Y),calc.r2,Y,pred)
    prior_covar = mmbr::create_mash_prior(sample_data =
                                              list(X=X,
                                                   Y=Y,
                                                   residual_variance =
                                                       cov(Y)),
                                          max_mixture_len=-1)
    m <- mmbr::msusie(X = X,
                      Y = Y,
                      prior_variance = prior_covar,
                      precompute_covariances = T,
                      L = L,
                      intercept = T,
                      standardize = T)

    modelList = list()
    for (i in 1:ncol(Y)){

        mod = tibble::tibble(SNP = colnames(X),
                             Weight = m$coef[-1,i])
        mod = subset(mod,Weight != 0)
        modelList = rlist::list.append(modelList,
                                       list(Transcript = colnames(Y)[i],
                                            Model = mod,
                                            R2 = r2.vec[i]))

    }

    return(modelList)


}
