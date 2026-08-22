recover_dropout_expression <- function(x, mask, membership, alpha=0.75, adjacency=NULL){
  out <- as.matrix(x)
  for(i in seq_len(nrow(x))){
    for(j in seq_len(ncol(x))){
      if(mask[i,j]){
        id <- which(membership==membership[j] & seq_along(membership)!=j)
        cell_pred <- mean(x[i,id],na.rm=TRUE)
        gene_pred <- if(is.null(adjacency)) cell_pred else sum(adjacency[i,]*x[,j])/sum(adjacency[i,])
        out[i,j] <- alpha*cell_pred+(1-alpha)*gene_pred
      }
    }
  }
  out
}

DropoutKiller <- function(x,membership,rank=20,threshold=0.95,alpha=0.75,adjacency=NULL){
  score <- local_alra_score(x,membership,rank)
  mask <- select_dropout_mask(score,threshold)
  list(expression=recover_dropout_expression(x,mask,membership,alpha,adjacency),score=score,mask=mask)
}

dropout_killer <- DropoutKiller
