local_alra_score <- function(x, membership, rank=20){
  stopifnot(inherits(x,"matrix") || inherits(x,"dgCMatrix"))
  score <- matrix(0,nrow(x),ncol(x),dimnames=dimnames(x))
  for(k in unique(membership)){
    id <- which(membership==k)
    y <- as.matrix(x[,id,drop=FALSE])
    fit <- irlba::irlba(y, nv=min(rank,ncol(y)-1), nu=min(rank,nrow(y)-1))
    yh <- fit$u %*% diag(fit$d) %*% t(fit$v)
    tau <- apply(yh,1,quantile,probs=0.05)
    z <- sweep(yh,1,tau)
    score[,id] <- 1/(1+exp(-z))
  }
  score
}

select_dropout_mask <- function(score, threshold=0.95) score>threshold
