.dk_atac_assay <- function(object, atac_assay=NULL) {
  assays<-names(object@assays)
  if(!is.null(atac_assay)) {
    if(length(atac_assay)!=1L||!atac_assay%in%assays)stop('atac_assay not found',call.=FALSE)
    return(atac_assay)
  }
  hits<-assays[vapply(object@assays,inherits,logical(1),what='ChromatinAssay') |
                tolower(assays)%in%c('atac','peaks')]
  if(length(hits)>1L)stop('multiple ATAC assays found; specify atac_assay',call.=FALSE)
  if(length(hits))hits else NULL
}

.dk_seurat_layer <- function(object,assay,layer='counts') {
  if(inherits(object[[assay]],'Assay5')) {
    available<-SeuratObject::Layers(object[[assay]],search=NA)
    if(!layer%in%available)stop('assay ',assay,' needs one joined ',layer,' layer; use JoinLayers first',call.=FALSE)
    SeuratObject::LayerData(object,assay=assay,layer=layer)
  } else if (utils::packageVersion('SeuratObject') >= '5.0.0') {
    SeuratObject::GetAssayData(object,assay=assay,layer=layer)
  } else SeuratObject::GetAssayData(object,assay=assay,slot=layer)
}

# Faithful kernel transformation from SuperCell 2.0 ComputeMultimodalKnn:
# affinity = 1 - 2 * weighted.nn distance^2; A + t(A); no self edges.
.dk_wnn_affinity <- function(index,distance,cells) {
  n<-length(cells)
  if(!identical(dim(index),dim(distance))||nrow(index)!=n||anyNA(index)||
     any(index<1L|index>n)||any(!is.finite(distance)))stop('invalid weighted.nn neighbors',call.=FALSE)
  a<-Matrix::sparseMatrix(i=rep(seq_len(n),each=ncol(index)),j=as.vector(t(index)),
      x=as.vector(t(1-2*distance^2)),dims=c(n,n),dimnames=list(cells,cells))
  a<-a+Matrix::t(a);diag(a)<-0
  if(any(a@x<0))a@x<-pmax(a@x,1e-16)
  methods::as(Matrix::drop0(a),'dgCMatrix')
}

#' Rebuild RNA PCA, ATAC LSI and SuperCell 2.0 WNN hierarchy
#' @param object Seurat object containing RNA and ATAC counts.
#' @param rna_assay RNA assay name.
#' @param atac_assay ATAC assay name or NULL for unambiguous automatic detection.
#' @param group Named broad cell-class labels. No fine annotation is used automatically.
#' @param gamma Cells per metacell.
#' @param npcs Maximum PCA and LSI dimensions; LSI component one is always excluded.
#' @param k_nn WNN neighbor count within each broad class.
#' @param seed Random seed.
#' @return Internal geometry with RNA embedding, WNN affinity, retained Walktrap tree and provenance.
#' @export
build_wnn_supercell <- function(object,rna_assay='RNA',atac_assay=NULL,group=NULL,
                                gamma=150,npcs=40L,k_nn=30L,seed=1L) {
  if(!requireNamespace('Seurat',quietly=TRUE)||!requireNamespace('Signac',quietly=TRUE))
    stop('Seurat and Signac are required to rebuild RNA+ATAC WNN',call.=FALSE)
  if(!inherits(object,'Seurat'))stop('object must be Seurat',call.=FALSE)
  atac_assay<-.dk_atac_assay(object,atac_assay)
  if(is.null(atac_assay)||!rna_assay%in%names(object@assays)||identical(rna_assay,atac_assay))
    stop('distinct RNA and ATAC assays are required',call.=FALSE)
  for(v in list(npcs,k_nn))if(length(v)!=1L||!is.finite(v)||v<2)stop('npcs and k_nn must be >= 2',call.=FALSE)
  rna<-.dk_seurat_layer(object,rna_assay);atac<-.dk_seurat_layer(object,atac_assay)
  cells<-colnames(rna)
  if(!setequal(cells,colnames(atac)))stop('RNA and ATAC must contain the same paired cells',call.=FALSE)
  if(!identical(cells,colnames(object)))object<-object[,cells]
  group<-.dk_align_vector(group,cells,'group')
  if(is.null(group))group<-rep('all',length(cells))
  if(anyNA(group))stop('broad groups cannot be missing',call.=FALSE)
  # Always recompute both reductions and WNN; old reductions/graphs are not read.
  SeuratObject::DefaultAssay(object)<-rna_assay
  object<-Seurat::NormalizeData(object,verbose=FALSE)
  object<-Seurat::FindVariableFeatures(object,nfeatures=min(2000L,nrow(rna)),verbose=FALSE)
  features<-SeuratObject::VariableFeatures(object)
  nr<-min(as.integer(npcs),length(cells)-1L,length(features)-1L)
  if(nr<2L)stop('insufficient RNA features/cells for PCA',call.=FALSE)
  object<-Seurat::ScaleData(object,features=features,verbose=FALSE)
  object<-Seurat::RunPCA(object,features=features,npcs=nr,reduction.name='dk_pca',seed.use=seed,verbose=FALSE)
  SeuratObject::DefaultAssay(object)<-atac_assay
  object<-Signac::RunTFIDF(object,verbose=FALSE)
  object<-Signac::FindTopFeatures(object,min.cutoff='q0',verbose=FALSE)
  na<-min(as.integer(npcs),length(cells)-1L,length(SeuratObject::VariableFeatures(object))-1L)
  if(na<2L)stop('ATAC needs at least two LSI dimensions; LSI1 is never used as a fallback',call.=FALSE)
  set.seed(seed)
  object<-Signac::RunSVD(object,n=na,reduction.name='dk_lsi',verbose=FALSE)
  pca<-SeuratObject::Embeddings(object[['dk_pca']])[cells,,drop=FALSE]
  nr<-ncol(pca);na<-ncol(SeuratObject::Embeddings(object[['dk_lsi']]))
  if(na<2L)stop('ATAC LSI has fewer than two components',call.=FALSE)
  all_i<-all_j<-integer();all_w<-numeric();small<-character()
  for(s in unique(group)) {
    ids<-which(group==s)
    if(length(ids)<4L) {small<-c(small,as.character(s));next}
    sub<-object[,cells[ids]];kk<-min(as.integer(k_nn),length(ids)-2L)
    set.seed(seed)
    sub<-Seurat::FindMultiModalNeighbors(sub,reduction.list=list('dk_pca','dk_lsi'),
      dims.list=list(seq_len(nr),seq.int(2L,na)),k.nn=kk,knn.range=min(200L,length(ids)-1L),
      weighted.nn.name='dk_weighted.nn',knn.graph.name='dk_wknn',snn.graph.name='dk_wsnn',verbose=FALSE)
    nn<-sub@neighbors[['dk_weighted.nn']]
    a<-.dk_wnn_affinity(nn@nn.idx,nn@nn.dist,cells[ids]);e<-Matrix::summary(a)
    all_i<-c(all_i,ids[e$i]);all_j<-c(all_j,ids[e$j]);all_w<-c(all_w,e$x)
  }
  if(length(small))warning('Broad strata with fewer than four cells retain singleton memberships and neutral neighbor support: ',paste(small,collapse=', '),call.=FALSE)
  a<-Matrix::sparseMatrix(i=all_i,j=all_j,x=all_w,dims=c(length(cells),length(cells)),dimnames=list(cells,cells))
  list(embedding=pca,affinity=a,membership_fit=.dk_lean_graph_membership(a,group,gamma),
       provenance=list(method='SuperCell2 ComputeMultimodalKnn kernel + Walktrap',
         reduction.list=c('PCA','LSI'),dims.list=list(seq_len(nr),seq.int(2L,na)),
         rebuilt=TRUE,rna_assay=rna_assay,atac_assay=atac_assay,k_nn=k_nn,
         small_strata=small,seed=seed))
}
