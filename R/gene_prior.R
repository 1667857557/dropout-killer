gene_prior_prediction <- function(x, adjacency, normalize=TRUE){
  if(normalize){
    s <- rowSums(adjacency)
    s[s==0] <- 1
    adjacency <- adjacency/s
  }
  as.matrix(adjacency%*%x)
}

combine_gene_prior <- function(x, ppi=NULL, pathway=NULL, grn=NULL, weights=c(0.33,0.33,0.34)){
  priors <- list(ppi,pathway,grn)
  priors <- priors[!vapply(priors,is.null,logical(1))]
  if(!length(priors)) return(x)
  Reduce(`+`,Map(function(a,w)w*gene_prior_prediction(x,a),priors,weights[seq_along(priors)]))
}
