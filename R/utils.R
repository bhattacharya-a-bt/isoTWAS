"%^%" <- function(x, n) {
    with(eigen(x), vectors %*% (values^n * t(vectors)))
    }
