selective_recovery <- function(x, mask, cell_prediction, gene_prediction, alpha=0.75){
  out <- x
  pred <- alpha*cell_prediction+(1-alpha)*gene_prediction
  out[mask] <- pred[mask]
  out
}

run_dropout_killer <- function(x, membership, embedding=NULL, adjacency=NULL, rank=20, threshold=0.95, alpha=0.75, sigma=1, k=50){
  score <- local_alra_score(x,membership,rank)
  mask <- select_dropout_mask(score,threshold)
  if(is.null(embedding)) embedding <- matrix(0,ncol(x),1)
  cell_pred <- weighted_neighbor_prediction(x,embedding,membership,sigma,k)
  gene_pred <- if(is.null(adjacency)) cell_pred else gene_prior_prediction(x,adjacency)
  list(expression=selective_recovery(x,mask,cell_pred,gene_pred,alpha),mask=mask,score=score,cell_prediction=cell_pred,gene_prediction=gene_pred)
}
