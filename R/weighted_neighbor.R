weighted_neighbor_prediction <- function(x, embedding, membership, sigma=1, k=50){
  stopifnot(ncol(x)==nrow(embedding),length(membership)==ncol(x))
  out <- matrix(0,nrow(x),ncol(x),dimnames=dimnames(x))
  for(i in seq_len(ncol(x))){
    idx <- which(membership==membership[i] & seq_len(ncol(x))!=i)
    if(length(idx)>k){
      d <- rowSums((embedding[idx,,drop=FALSE]-embedding[i,])^2)
      idx <- idx[order(d)[seq_len(k)]]
      d <- d[order(d)[seq_len(k)]]
    }else d <- rowSums((embedding[idx,,drop=FALSE]-embedding[i,])^2)
    if(length(idx)){
      w <- exp(-d/(sigma^2));w <- w/sum(w)
      out[,i] <- as.vector(as.matrix(x[,idx,drop=FALSE])%*%w)
    }
  }
  out
}
