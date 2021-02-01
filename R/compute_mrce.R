#' Multivariate regression with covariance estimate
#'
#' The function trains multivariate lasso with a given covariance estimate
#'
#' @param X matrix, design matrix of SNP dosages
#' @param Y matrix, matrix of G isoform expression across columns
#' @param lambda vector, lambda penalty vector for LASSO to tune on
#' @param Omega matrix, precision matrix of Y
#' @param nfolds int, number of CV folds
#' @param tol.in numeric, tolerance for objective difference
#' @param maxit.in int, maximum number of iteractions
#' @param verbose logical
#' @param seed int, random seed
#'
#' @return CV MRCE fit
#'
#' @export
compute_mrce = function(X,
                        Y,
                        lambda = NULL,
                        nlambda = 50,
                        Omega,
                        nfolds = 5,
                        tol.in,
                        maxit.in,
                        verbose,
                        seed){

    if (is.null(lambda)){
        lambda = 10^(seq(2,-30,length.out = nlambda))
        }

    set.seed(seed)
    train.folds = caret::createFolds(1:nrow(Y),
                                     k = nfolds,
                                     returnTrain = T)

    get_cvm <- function(l,
                        nfolds){

        pred = matrix(nrow = nrow(Y),
                      ncol = ncol(Y))
        for (tr in 1:nfolds){

            Y.tr = Y[train.folds[[tr]],]
            X.tr = X[train.folds[[tr]],]
            X.test = X[-train.folds[[tr]],]

            cur.B = compute_fixed(X.tr,
                                  Y.tr,
                                  l,
                                  Omega,
                                  tol.in,
                                  maxit.in,
                                  silent = T)$Bhat

            pred[-train.folds[[tr]],] = X.test %*% cur.B

        }

        return(sapply(1:ncol(pred),calc.r2,Y,pred))

    }

    cvm = pbapply::pbsapply(lambda,get_cvm,nfolds=nfolds)

    final.model = compute_fixed(X = X,
                                Y = Y,
                                l = lambda[which.max(colMeans(cvm))],
                                Omega = Omega,
                                tol.in = tol.in,
                                maxit.in = maxit.in,
                                silent = verbose)


    modelList = list()
    for (i in 1:ncol(Y)){

        mod = tibble::tibble(SNP = colnames(X),
                             Weight = final.model$Bhat[,i])
        mod = subset(mod,Weight != 0)
        modelList = rlist::list.append(modelList,
                                       list(Transcript = colnames(Y)[i],
                                            Model = mod,
                                            R2 = cvm[i,
                                                     which.max(colMeans(cvm))]))

    }

    return(modelList)




}
