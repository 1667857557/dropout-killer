#' Calibrate Lean hierarchy detection using artificial UMI thinning
#' @param counts Raw integer RNA counts (genes by cells).
#' @param group Broad cell classes used to define proxy negatives and hard boundaries.
#' @param q UMI retention probability, not the fraction of matrix entries masked.
#' @param seeds Independent thinning seeds. All resulting events are used.
#' @param fpr Training proxy-FPR target.
#' @param negative_max Maximum prevalence in the negative's own broad class.
#' @param negative_other_min Minimum prevalence in another broad class.
#' @param gamma Cells per SuperCell membership.
#' @param k_knn RNA graph neighbor count.
#' @param neighbor_k Maximum hierarchy support neighbors.
#' @param rank ALRA rank or auto.
#' @param geometry_builder Optional function(counts) returning freshly rebuilt WNN geometry.
#' @return DropoutKillerLeanModel; use separate masks/data for evaluation.
#' @export
calibrate_lean_detector <- function(counts, group, q=.5, seeds=10001:10003,
                                    fpr=.01, negative_max=.005, negative_other_min=.2,
                                    gamma=150, k_knn=5L, neighbor_k=30L, rank='auto',
                                    geometry_builder=NULL) {
  counts<-.dk_validate_expression(counts);nm<-.dk_names(counts);dimnames(counts)<-list(nm$genes,nm$cells)
  counts<-methods::as(Matrix::Matrix(counts,sparse=TRUE),'dgCMatrix')
  if(any(abs(counts@x-round(counts@x))>1e-8))stop('calibration requires raw integer UMI counts',call.=FALSE)
  group<-.dk_align_vector(group,nm$cells,'group',allow_null=FALSE)
  if(anyNA(group)||length(unique(group))<2L)stop('calibration needs at least two broad classes for proxy negatives; supply a fitted lean_model otherwise',call.=FALSE)
  if(length(q)!=1L||!is.finite(q)||q<=0||q>=1)stop('q must be in (0,1)',call.=FALSE)
  if(!length(seeds)||anyNA(seeds)||any(!is.finite(seeds)))stop('provide finite calibration seeds',call.=FALSE)
  for(v in list(negative_max,negative_other_min))if(length(v)!=1L||!is.finite(v)||v<0||v>1)stop('negative prevalence bounds must be in [0,1]',call.=FALSE)
  lev<-unique(group)
  prevalence<-vapply(lev,function(s)Matrix::rowMeans(counts[,group==s,drop=FALSE]>0),numeric(nrow(counts)))
  negatives<-list()
  for(s in seq_along(lev)) {
    good<-which(prevalence[,s]<=negative_max & apply(prevalence[,-s,drop=FALSE],1,max)>=negative_other_min)
    for(g in good) {
      cc<-which(group==lev[s] & as.numeric(counts[g,])==0)
      if(length(cc))negatives[[length(negatives)+1L]]<-cbind(g,cc)
    }
  }
  negative<-do.call(rbind,negatives)
  if(is.null(negative))stop('no valid proxy negatives: supply an independently fitted lean_model or explicitly revise negative_max; no silent relaxation',call.=FALSE)
  features<-truth<-vector('list',length(seeds));methods_seen<-character()
  for(j in seq_along(seeds)) {
    set.seed(seeds[j]);masked<-counts
    masked@x<-as.numeric(stats::rbinom(length(counts@x),counts@x,q))
    lost<-which(masked@x==0 & counts@x>0)
    cc<-rep.int(seq_len(ncol(counts)),diff(counts@p))
    positive<-cbind(counts@i[lost]+1L,cc[lost])
    if(!nrow(positive))stop('thinning produced no positive events',call.=FALSE)
    masked<-Matrix::drop0(masked)
    # Empty libraries are valid after thinning; keep them at zero on this scale.
    y<-masked%*%Matrix::Diagonal(x=1e4/pmax(Matrix::colSums(masked),1));y@x<-log1p(y@x)
    dimnames(y)<-dimnames(masked)
    sv<-.dk_lean_svd(y,rank,seeds[j]+300000L)
    geometry<-if(is.null(geometry_builder)).dk_lean_rna_geometry(sv$embedding,group,gamma,k_knn)else geometry_builder(masked)
    methods_seen<-c(methods_seen,.dk_lean_method(!is.null(geometry$affinity)))
    .dk_validate_lean_geometry(geometry,colnames(counts),group)
    nb<-.dk_lean_neighbors(geometry,neighbor_k)
    features[[j]]<-.dk_lean_features(y,rbind(positive,negative),geometry,sv,nb)
    truth[[j]]<-c(rep(TRUE,nrow(positive)),rep(FALSE,nrow(negative)))
  }
  if(length(unique(methods_seen))!=1L)stop('calibration geometry changed modality',call.=FALSE)
  model<-fit_lean_detector(do.call(rbind,features),unlist(truth),fpr,methods_seen[1])
  model$calibration<-list(q=q,seeds=seeds,negative_max=negative_max,negative_other_min=negative_other_min,
                          geometry='rebuilt from each thinned RNA matrix',all_events=TRUE)
  model
}
