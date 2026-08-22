recover_dropout_expression <- function(x, mask, membership, embedding=NULL, alpha=0.75, adjacency=NULL, k=50, sigma=1){
  out <- as.matrix(x)
  cell_pred <- if(is.null(embedding)) {
    sapply(seq_len(ncol(x)), function(j) rowMeans(x[,which(membership==membership[j] & seq_len(ncol(x))!=j),drop=FALSE]))
  } else weighted_neighbor_prediction(x,embedding,membership,k=k,sigma=sigma)
  if(is.vector(cell_pred)) cell_pred<-matrix(cell_pred,ncol=1)
  gene_pred <- if(is.null(adjacency)) cell_pred else gene_prior_prediction(x,adjacency)
  pred <- combine_recovery_prediction(cell_pred,gene_pred,alpha)
  out[mask] <- pred[mask]
  out
}

DropoutKiller <- function(x,membership,embedding=NULL,rank=20,threshold=0.95,alpha=0.75,adjacency=NULL,k=50,sigma=1){
  score <- local_alra_score(x,membership,rank)
  mask <- select_dropout_mask(score,threshold)
  list(expression=recover_dropout_expression(x,mask,membership,embedding,alpha,adjacency,k,sigma),score=score,mask=mask)
}

dropout_killer <- DropoutKiller
