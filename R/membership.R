dropout_membership <- function(object, reduction="pca", k=100, seed=1){
  if(!requireNamespace("SuperCell", quietly=TRUE)) stop("SuperCell required")
  set.seed(seed)
  z <- Seurat::Embeddings(object, reduction=reduction)
  mem <- SuperCell::SCimplify(z, k=k, seed=seed)
  mem
}
