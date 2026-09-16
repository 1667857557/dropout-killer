lean_fixture <- function() {
  set.seed(901)
  x<-matrix(rpois(24*24,2),24,24,dimnames=list(paste0('g',1:24),paste0('c',1:24)))
  z<-matrix(rnorm(24*3),24,3,dimnames=list(colnames(x),NULL))
  group<-rep(c('B','T'),each=12)
  list(x=x,z=z,group=group)
}
lean_test_model <- function(wnn=FALSE) {
  set.seed(902)
  f<-data.frame(z=rnorm(600),local_hierarchy=rnorm(600),H=rnorm(600))
  truth<-runif(600)<plogis(f$z+.5*f$local_hierarchy+.25*f$H)
  fit_lean_detector(f,truth,method=.dk_lean_method(wnn))
}

test_that('WNN affinity follows the actual SuperCell 2.0 kernel', {
  idx<-matrix(c(2,3,1,3,1,2),3,2,byrow=TRUE)
  d<-matrix(c(.2,.4,.3,.5,.4,.6),3,2,byrow=TRUE)
  a<-.dk_wnn_affinity(idx,d,c('a','b','c'))
  expected<-matrix(0,3,3)
  for(i in 1:3)for(j in 1:2)expected[i,idx[i,j]]<-1-2*d[i,j]^2
  expect_equal(unname(as.matrix(a)),expected+t(expected))
  expect_equal(diag(as.matrix(a)),c(a=0,b=0,c=0))
})

test_that('Walktrap keeps every cell and respects broad classes and components', {
  d<-lean_fixture();geo<-.dk_lean_rna_geometry(d$z,d$group,gamma=4,k_knn=3)
  fit<-geo$membership_fit
  expect_true(all(vapply(split(d$group,fit$membership),function(x)length(unique(x))==1,logical(1))))
  expect_true(all(fit$strata$target_memberships==3))
  a<-Matrix::sparseMatrix(i=c(1,2),j=c(2,1),x=1,dims=c(4,4),dimnames=list(letters[1:4],letters[1:4]))
  f<-.dk_lean_graph_membership(a,rep('B',4),gamma=150)
  expect_equal(length(f$membership),4)
  expect_equal(length(unique(f$membership)),3)
})

test_that('hierarchy support excludes self, uses own membership first and borrows on shortage', {
  d<-lean_fixture();geo<-.dk_lean_rna_geometry(d$z,d$group,gamma=4,k_knn=4)
  nb<-.dk_lean_neighbors(geo,k=5);e<-Matrix::summary(nb$W);m<-geo$membership_fit$membership
  expect_true(all(e$i!=e$j));expect_true(all(d$group[e$i]==d$group[e$j]))
  expect_equal(as.numeric(Matrix::rowSums(nb$W)),rep(1,24),tolerance=1e-12)
  for(i in 1:24) {
    own<-which(m==m[i] & seq_along(m)!=i)
    donors<-e$j[e$i==i]
    if(length(own)>=5)expect_true(all(donors%in%own)) else expect_true(all(own%in%donors))
  }
  expect_gt(sum(nb$borrowed),0)
})

test_that('Lean features equal scalar benchmark definitions', {
  d<-lean_fixture();x<-.dk_alra_library_log(d$x);sv<-.dk_lean_svd(x,rank=3,seed=4)
  geo<-.dk_lean_rna_geometry(d$z,d$group,4,4);nb<-.dk_lean_neighbors(geo,5)
  ij<-which(x==0,arr.ind=TRUE);f<-.dk_lean_features(x,ij,geo,sv,nb)
  m<-geo$membership_fit$membership
  for(a in seq_len(nrow(ij))) {
    g<-ij[a,1];c<-ij[a,2];own<-which(m==m[c]);grp<-which(d$group==d$group[c])
    H<-length(own)/(length(own)+50)*(qlogis((sum(x[g,own]>0)+.5)/(length(own)+1))-qlogis((sum(x[g,grp]>0)+.5)/(length(grp)+1)))
    local<-qlogis((sum(nb$W[c,]*as.numeric(x[g,]>0))+.01)/1.02)
    expect_equal(f$H[a],H,tolerance=1e-12)
    expect_equal(f$local_hierarchy[a],local,tolerance=1e-12)
    lr<-as.numeric(sv$U%*%sv$V[g,]);sig<-sd(lr[lr<0])
    if(!is.finite(sig)||sig<=0)sig<-sd(lr)
    if(!is.finite(sig)||sig<1e-8)sig<-1e-8
    expect_equal(f$z[a],(lr[c]-abs(quantile(lr,.001,names=FALSE)))/sig,tolerance=1e-12)
  }
})

test_that('model fitting and strict proxy-FPR threshold preserve ties', {
  model<-lean_test_model();expect_s3_class(model,'DropoutKillerLeanModel')
  set.seed(902);f<-data.frame(z=rnorm(600),local_hierarchy=rnorm(600),H=rnorm(600))
  truth<-runif(600)<plogis(f$z+.5*f$local_hierarchy+.25*f$H)
  scores<-as.numeric(model$beta[1]+as.matrix(f)%*%model$beta[-1])
  expect_lte(mean(scores[!truth]>model$threshold),.01)
  expect_error(.dk_validate_lean_model(model,.dk_lean_method(TRUE)),'RNA/WNN')
  f$H<-0;expect_true(fit_lean_detector(f,truth)$converged)
})

test_that('new default calls only zeros and leaves observed recovery values unchanged', {
  d<-lean_fixture();m<-lean_test_model();m$threshold<--100
  res<-dropout_killer(d$x,d$z,group=d$group,lean_model=m,gamma=4,k_knn=4,
                     rank=3,recovery_method='neighbor',neighbor_k=5)
  expect_identical(res$settings$detection_method,.dk_lean_method())
  expect_true(validate_dropout_result(res,d$x)$valid)
  expect_true(all(d$x[cbind(res$events$i,res$events$j)]==0))
  m$threshold<-1e10
  res<-dropout_killer(d$x,d$z,group=d$group,lean_model=m,gamma=4,rank=3,recovery_method='neighbor')
  expect_equal(nrow(res$events),0)
  expect_error(dropout_killer(d$x,d$z,group=d$group,membership=rep(1,24),lean_model=m),'retained hierarchy')
  m$threshold<--100
  p1<-dropout_killer(d$x,d$z,group=d$group,lean_model=m,gamma=4,k_knn=4,rank=3,
                     factor_features=12,min_feature_observed=2,min_target_observed=2,
                     factor_crossfit_folds=2)
  expect_true(validate_dropout_result(p1,d$x)$valid)
  expect_identical(p1$settings$recovery_method,'p1_stabilized_state')
})

test_that('calibration fails clearly without proxy negatives and rebuilds geometry', {
  d<-lean_fixture()
  expect_error(calibrate_lean_detector(d$x,rep('B',24)),'two broad classes')
  expect_error(calibrate_lean_detector(d$x,d$group),'no valid proxy negatives')
  d$x[1:4,1:12]<-0;d$x[5:8,13:24]<-0
  calls<-0L
  builder<-function(counts) {
    calls<<-calls+1L
    expect_lt(sum(counts),sum(d$x))
    y<-counts%*%Matrix::Diagonal(x=1e4/pmax(Matrix::colSums(counts),1));y@x<-log1p(y@x)
    dimnames(y)<-dimnames(counts)
    .dk_lean_rna_geometry(.dk_lean_svd(y,3,42)$embedding,d$group,4,4)
  }
  model<-calibrate_lean_detector(d$x,d$group,seeds=c(31,32),rank=3,gamma=4,geometry_builder=builder)
  expect_equal(calls,2L);expect_true(model$converged);expect_true(model$calibration$all_events)
})

test_that('Seurat RNA default does not require stored reductions', {
  skip_if_not_installed('Seurat')
  d<-lean_fixture();obj<-SeuratObject::CreateSeuratObject(d$x)
  obj$broad<-d$group
  out<-dropout_killer_seurat(obj,group_by='broad',lean_model=lean_test_model(),return_result=TRUE,
                            gamma=4,rank=3,recovery_method='neighbor')
  expect_identical(out$result$settings$detection_method,.dk_lean_method())
  expect_true(validate_dropout_result(out$result,d$x,tolerance=1e-12)$valid)
  expect_error(dropout_killer_seurat(obj,modality='rna',detection_method=.dk_lean_method(TRUE),lean_model=lean_test_model()),'conflicts')
})

test_that('paired multiome rebuilds PCA and LSI2+ and automatically uses WNN', {
  skip_if_not_installed('Seurat');skip_if_not_installed('Signac')
  set.seed(910);cells<-paste0('c',1:60)
  rna<-matrix(rpois(120*60,2),120,60,dimnames=list(paste0('g',1:120),cells))
  atac<-matrix(rpois(160*60,.7),160,60,dimnames=list(paste0('p',1:160),cells))
  obj<-SeuratObject::CreateSeuratObject(rna)
  obj[['ATAC']]<-SeuratObject::CreateAssayObject(counts=atac)
  obj$broad<-rep(c('B','T'),each=30)
  SeuratObject::DefaultAssay(obj)<-'ATAC'
  geo<-build_wnn_supercell(obj,group=obj$broad,gamma=5,npcs=5,k_nn=5)
  expect_identical(geo$provenance$dims.list[[2]],2:5)
  expect_true(geo$provenance$rebuilt)
  e<-Matrix::summary(geo$affinity)
  expect_true(all(obj$broad[e$i]==obj$broad[e$j]));expect_true(all(e$i!=e$j))
  expect_equal(geo$affinity,Matrix::t(geo$affinity),ignore_attr=TRUE)
  nb<-.dk_lean_neighbors(geo,k=3);ne<-Matrix::summary(nb$W)
  expect_true(all(geo$affinity[cbind(ne$i,ne$j)]>0))
  mem<-geo$membership_fit$membership
  for(i in seq_along(cells)) {
    own<-which(mem==mem[i] & as.numeric(geo$affinity[i,])>0)
    if(length(own)>=3)expect_true(all(ne$j[ne$i==i]%in%own))
  }
  perm<-geo;perm$embedding<-perm$embedding[rev(seq_along(cells)),,drop=FALSE]
  expect_error(.dk_validate_lean_geometry(perm,cells,obj$broad),'alignment')
  m<-lean_test_model(TRUE);m$threshold<-1e10
  out<-dropout_killer_seurat(obj,group_by='broad',lean_model=m,return_result=TRUE,wnn_npcs=5,wnn_k=5,
    gamma=5,rank=3,recovery_method='neighbor')
  expect_identical(out$result$settings$detection_method,.dk_lean_method(TRUE))
  expect_identical(out$result$settings$wnn$rna_assay,'RNA')
  expect_identical(out$result$settings$wnn$dims.list[[2]],2:5)
})
