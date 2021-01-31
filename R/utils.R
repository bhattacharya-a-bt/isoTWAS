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
