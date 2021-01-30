### SOME SHARED EFFECTS
set.seed(1218)
N = 150
R = 50
P = 700
G = 5

Posdef <- function (n,
                    ev = runif(n, 0, 10)) {
    Z <- matrix(ncol=n, rnorm(n^2))
    decomp <- qr(Z)
    Q <- qr.Q(decomp)
    R <- qr.R(decomp)
    d <- diag(R)
    ph <- d / abs(d)
    O <- Q %*% diag(ph)
    Z <- t(O) %*% diag(ev) %*% O
    return(Z)
}

maf = runif(P, min = 0.05, max = 0.5)
snps = sapply(maf,function(p) rbinom(N,size = 2,prob = p))
beta = MASS::mvrnorm(n = P,
                     mu = rep(0,G),
                     Sigma = Posdef(G))
causals_shared = rbinom(P,1,.05)
causals_all = replicate(5,causals_shared + rbinom(P,1,.05))
causals_all[causals_all > 1] = 1
beta[causals_all == 0] = 0
Y = snps %*% beta + matrix(rnorm(N*G),nrow=N)

Y.list = list()
for (i in 1:R){

    Y.list = rlist::list.append(Y.list,
                                jitter(Y,amount=5))

}

colnames(Y) = paste0('Isoform',1:G)
colnames(snps) = paste0('SNP',1:P)

Y.full = abind::abind(Y.list,along = 1)
total = as.data.frame(cbind(Y.full,rep(1:N,R)))
colnames(total) = c(paste0('Isoform',1:5),
                    'Sample')
total$Sample = as.factor(total$Sample)

Omega = huge::huge(Y,method='glasso')$icov[[10]]

