\name{build_wnn_supercell}
\alias{build_wnn_supercell}
\title{Rebuild RNA PCA, ATAC LSI and SuperCell 2.0 WNN hierarchy}
\description{Rebuild RNA PCA, ATAC LSI and SuperCell 2.0 WNN hierarchy.}
\usage{
build_wnn_supercell(object,rna_assay='RNA',atac_assay=NULL,group=NULL,
                                gamma=150,npcs=40L,k_nn=30L,seed=1L)
}
\arguments{
\item{object}{Seurat object containing RNA and ATAC counts.}
\item{rna_assay}{RNA assay name.}
\item{atac_assay}{ATAC assay name or NULL for unambiguous automatic detection.}
\item{group}{Named broad cell-class labels. No fine annotation is used automatically.}
\item{gamma}{Cells per metacell.}
\item{npcs}{Maximum PCA and LSI dimensions; LSI component one is always excluded.}
\item{k_nn}{WNN neighbor count within each broad class.}
\item{seed}{Random seed.}
}
\value{Internal geometry with RNA embedding, WNN affinity, retained Walktrap tree and provenance.}
