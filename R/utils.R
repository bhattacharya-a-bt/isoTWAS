"%^%" <- function(x, n) {
    with(eigen(x), vectors %*% (values^n * t(vectors)))
    }

.onUnload <- function (libpath) { library.dynam.unload("isoTWAS", libpath)}

calc.mse <- function(obs, pred){
    if(is.vector(obs)) obs <- as.matrix(obs)
    if(is.vector(pred)) pred <- as.matrix(pred)

    n <- nrow(obs)
    rss <- colSums((obs - pred)^2, na.rm = TRUE)
    rss/n
}

calc.r2 <- function(i,obs,pred){

    summary(lm(obs[,i] ~ pred[,i]))$adj.r.sq


}


get_best <- function(list_mods,G=G){
    r2 = sapply(list_mods, function(y) sapply(y,function(x) x$R2))
    bundle = cbind(apply(r2,1,which.max),1:G)
    return(lapply(1:G,
                  function(i) list_mods[[bundle[i,1]]][[bundle[i,2]]]))
}
