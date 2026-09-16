lean_fixture <- function() {
  set.seed(901)
  x<-matrix(rpois(24*24,2),24,24,dimnames=list(paste0('g',1:24),paste0('c',1:24)))
  z<-matrix(rnorm(24*3),24,3,dimnames=list(colnames(x),paste0('PC_',1:3)))
  list(x=x,z=z,group=rep(c('B','T'),each=12))
}
lean_test_model <- function() {
  set.seed(902)
  f<-data.frame(z=rnorm(600),local_hierarchy=rnorm(600),H=rnorm(600))
  truth<-runif(600)<plogis(f$z+.5*f$local_hierarchy+.25*f$H)
  fit_lean_detector(f,truth)
}

test_that('RNA opt-out and historical routes ignore ambiguous ATAC assays', {
  skip_if_not_installed('Seurat')
  d<-lean_fixture()
  obj<-SeuratObject::CreateSeuratObject(Matrix::Matrix(d$x,sparse=TRUE))
  obj[['ATAC']]<-SeuratObject::CreateAssayObject(counts=Matrix::Matrix(d$x,sparse=TRUE))
  obj[['peaks']]<-SeuratObject::CreateAssayObject(counts=Matrix::Matrix(d$x,sparse=TRUE))
  obj[['pca']]<-SeuratObject::CreateDimReducObject(embeddings=d$z,key='PC_',assay='RNA')
  obj$broad<-d$group
  rna<-suppressWarnings(dropout_killer_seurat(obj,modality='rna',atac_assay='unused',
    group_by='broad',lean_model=lean_test_model(),rank=3,gamma=4,
    recovery_method='neighbor',return_result=TRUE))
  expect_identical(rna$result$settings$detection_method,.dk_lean_method())
  for(method in c('alra_global_by_group','eb_zero_null','alra_quantile')) {
    old<-suppressWarnings(dropout_killer_seurat(obj,detection_method=method,
      group_by='broad',rank=3,gamma=4,min_cells=2,recovery_method='neighbor',return_result=TRUE))
    expect_identical(old$result$settings$detection_method,method)
  }
  expect_error(dropout_killer_seurat(obj,modality='auto'),'multiple ATAC assays')
  expect_error(dropout_killer_seurat(obj,modality='wnn'),'multiple ATAC assays')
})

test_that('automatic RNA calibration and inference use the same geometry pipeline', {
  d<-lean_fixture();d$x[1:4,1:12]<-0;d$x[5:8,13:24]<-0
  run<-function(z,model=NULL) dropout_killer(d$x,z,group=d$group,rank=3,gamma=4,
    k_knn=4,neighbor_k=5,seed=42,lean_model=model,
    lean_control=list(seeds=c(31,32)),recovery_method='neighbor')
  first<-run(d$z)
  # Deliberately change the external recovery embedding, preserving cell labels.
  other<-d$z[rev(seq_len(nrow(d$z))),,drop=FALSE];rownames(other)<-rownames(d$z)
  second<-run(other)
  expected<-.dk_lean_rna_geometry(.dk_lean_svd(.dk_alra_library_log(d$x),3,42)$embedding,d$group,4,4)
  expect_equal(unname(first$membership),unname(expected$membership_fit$membership))
  expect_equal(first$membership_fit$membership,second$membership_fit$membership)
  expect_equal(first$membership_fit$tree_indices,second$membership_fit$tree_indices)
  for(s in names(first$membership_fit$graphs))
    expect_equal(igraph::as_adjacency_matrix(first$membership_fit$graphs[[s]]),
                 igraph::as_adjacency_matrix(second$membership_fit$graphs[[s]]))
  expect_equal(first$detection$model,second$detection$model)
  expect_equal(first$detection$events,second$detection$events)
  expect_equal(first$mask,second$mask)
  reused<-run(other,first$detection$model)
  expect_equal(first$detection$events,reused$detection$events)
  expect_error(dropout_killer(d$x,d$z,group=d$group,lean_geometry=expected),
               'custom geometry requires a fitted lean_model')
  expect_error(dropout_killer(d$x,d$z,group=d$group,
    lean_control=list(geometry_builder=function(x)expected)),
    'custom geometry requires a fitted lean_model')
})
