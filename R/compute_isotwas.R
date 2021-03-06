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
#' @param omega_est character, 'replicates' or 'mean' to use Y.rep or Y
#' @param omega_nlambda int, number of omegas to generate
#' @param method character, vector of methods to use
#' @param predict_nlambda int, number of lambdas in MRCE
#' @param family character, glmnet family
#' @param scale logical, T/F to scale Y by Omega
#' @param alpha numeric, elastic net mixing parameter
#' @param nfolds int, number of CV folds
#' @param verbose logical
#' @param par logical, uses mclapply to parallelize model fit
#' @param n.cores int, number of parallel cores
#' @param tx_names vector, character vector of tx names - order of columns of Y
#' @param seed int, random seed
#' @param return_all logical, return R2 for all models?
#' @param tol.in numeric, tolerance for objective difference
#' @param maxit.in int, maximum number of iteractions
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
                            omega_est = 'replicates',
                            omega_nlambda = 10,
                            method = c('mrce_lasso','curds_whey',
                                       'multi_enet','mrce_enet',
                                       'mmbr','finemap',
                                       'univariate'),
                            predict_nlambda = 50,
                            family = 'gaussian',
                            scale = F,
                            alpha = 0.5,
                            nfolds = 5,
                            verbose = F,
                            par = F,
                            n.cores = NULL,
                            tx_names = NULL,
                            seed = NULL,
                            run_all = T,
                            return_all = F,
                            tol.in = 1e-3,
                            maxit.in = 1e3,
                            coverage = .9){

    ### CHECKS
    if (nrow(X) != nrow(Y)){
        stop('No. of rows of X =/= no. of rows of Y.')
    }

    N = nrow(Y)
    P = ncol(X)
    G = ncol(Y)

    if (is.null(R)){
        R = nrow(Y.rep)/nrow(Y)
    } else {
        if (nrow(Y.rep) != R*N){
            stop('No. of rows of Y.full =/= R*N')
        }
    }

    ### COMPUTE OMEGA
    if (verbose){print('Computing Omega')}
    omega_list = compute_omega(Y,
                               Y.rep,
                               R,
                               id,
                               method = omega_est,
                               nlambda = omega_nlambda,
                               verbose = verbose)

    if (run_all){
        print('run_all is set to T and all methods will be run.')
        method = c('mrce_lasso','curds_whey',
                   'multi_enet','mrce_enet',
                   'finemap','univariate')
        }

    if (is.null(seed)){
        seed = sample(1:100000,1)
    }

    all_models = list()

    if ('mrce_lasso' %in% method){
        if (verbose){print('Running mrce_lasso')}
        mrce_lasso = list()
        for (i in c(ceiling(omega_nlambda)/2,omega_nlambda)){
            mrce_lasso = rlist::list.append(mrce_lasso,
                                            compute_mrce(X = X,
                                                         Y = Y,
                                                         lambda = NULL,
                                                         nlambda =
                                                             predict_nlambda,
                                                         Omega =
                                                           omega_list$icov[[i]],
                                                         nfolds = nfolds,
                                                         tol.in = tol.in,
                                                         maxit.in = maxit.in,
                                                         verbose = verbose,
                                                         seed = seed))
        }
        all_models = rlist::list.append(all_models,get_best(mrce_lasso,
                                                            G = G))
    }

    if ('curds_whey' %in% method){
        if (verbose){print('Running curds_whey')}
        best_curds_whey = compute_curds_whey(X,
                                        Y,
                                        family = family,
                                        alpha = alpha,
                                        nfolds = nfolds,
                                        verbose = verbose,
                                        par = F,
                                        n.cores = NULL,
                                        tx_names = NULL,
                                        seed)
        all_models = rlist::list.append(all_models,best_curds_whey)
    }

    if ('multi_enet' %in% method){
        if (verbose){print('Running multi_enet')}
        best_multi_enet =
            multivariate_elasticnet(X = X,
                                    Y = Y,
                                    Omega =
                                        omega_list$icov[[omega_nlambda]],
                                    scale = scale,
                                    alpha = alpha,
                                    nfolds = nfolds,
                                    verbose = verbose,
                                    par = par,
                                    n.cores = n.cores,
                                    tx_names = tx_names,
                                    seed = seed)
        all_models = rlist::list.append(all_models,best_multi_enet)
    }

    if ('mmbr' %in% method){
        if (verbose){print('Running mmbr')}
        mmbr_mod = multivariate_mmbr(X = X,
                                     Y = Y,
                                     nfolds = nfolds,
                                     verbose = verbose,
                                     tx_names = tx_names,
                                     seed = seed)
        all_models = rlist::list.append(all_models,mmbr_mod)
    }

    if ('mrce_enet' %in% method){
        if (verbose){print('Running mrce_enet')}
        if (!verbose){
            sink(tempfile())
        }

        best_mrce_enet = compute_spring(X = X,
                                        Y = Y,
                                        Omega =
                                omega_list$icov[[length(omega_list$icov)]],
                                        nfolds = nfolds,
                                        tol.in = tol.in,
                                        maxit.in = maxit.in/10,
                                        verbose = verbose,
                                        seed = seed,
                                        par = par,
                                        n.cores = n.cores)

        if (!verbose){
            sink()
        }

        all_models = rlist::list.append(all_models,best_mrce_enet)
    }

    if ('finemap' %in% method){
        if (verbose){print('Running finemap')}
        best_finemap = compute_finemap_regress(X = X,
                                               Y = Y,
                                               Y.rep = Y.rep,
                                               R = R,
                                               id = id,
                                               nfolds = nfolds,
                                               verbose = verbose,
                                               tx_names = tx_names,
                                               coverage = coverage,
                                               seed = seed)
        all_models = rlist::list.append(all_models,best_finemap)
    }

    if ('univariate' %in% method){
        if (verbose){print('Running univariate')}
        uni_enet = univariate_elasticnet(X = X,
                                         Y = Y,
                                         Omega = omega_list[[omega_nlambda]],
                                         family = family,
                                         scale = scale,
                                         alpha = alpha,
                                         nfolds = nfolds,
                                         verbose = verbose,
                                         par = par,
                                         n.cores = n.cores,
                                         tx_names = tx_names,
                                         seed = seed)

        uni_blup = univariate_blup(X = X,
                                   Y = Y,
                                   Omega = omega_list[[omega_nlambda]],
                                   scale = scale,
                                   alpha = alpha,
                                   nfolds = nfolds,
                                   verbose = verbose,
                                   par = par,
                                   n.cores = n.cores,
                                   tx_names = tx_names,
                                   seed = seed)

        uni_susie = univariate_susie(X = X,
                                     Y = Y,
                                     Omega = omega_list[[omega_nlambda]],
                                     scale = scale,
                                     alpha = alpha,
                                     nfolds = nfolds,
                                     verbose = verbose,
                                     par = par,
                                     n.cores = n.cores,
                                     tx_names = tx_names,
                                     seed = seed)

        univariate = list(uni_enet,
                          uni_blup,
                          uni_susie)
        all_models = rlist::list.append(all_models,get_best(univariate,
                                                            G = G))

    }

    isotwas_mod = get_best(all_models,
                           G = G)

    if (return_all){
        r2 = sapply(all_models, function(y) sapply(y,function(x) x$R2))
        r2.df = as.data.frame(cbind(colnames(Y),r2))
        colnames(r2.df) = c('Transcript',method)
        isotwas_mod = list(Model = get_best(all_models,
                                            G = G),
                           R2 = r2.df)
        }

    return(isotwas_mod)


}
