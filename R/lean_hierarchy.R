# Lean features retain the benchmark definitions. A fitted calibration model is
# required: the three features are not combined with arbitrary fixed weights.
.dk_lean_method <- function(wnn = FALSE) {
  paste0('Supercell_hierarchy_Lean_membership', if (wnn) '_WNN' else '')
}

.dk_lean_svd <- function(x, rank = 'auto', seed = 1L) {
  A <- Matrix::t(x)
  kmax <- min(100L, min(dim(A)) - 1L)
  if (kmax < 1L) stop('Lean detection needs at least two genes and two cells', call. = FALSE)
  set.seed(seed)
  if (kmax >= min(dim(A)) - 1L) {
    f <- svd(as.matrix(A), nu = kmax, nv = kmax)
    f$d <- head(f$d, kmax)
  } else {
    f <- irlba::irlba(A, nv = kmax, nu = kmax,
                     work = max(2L * kmax + 10L, 160L), maxit = 2000)
  }
  if (identical(rank, 'auto')) {
    gaps <- -diff(f$d)
    tail <- if (length(gaps)) gaps[seq.int(min(80L, length(gaps)), length(gaps))] else numeric()
    s <- stats::sd(tail)
    hits <- if (is.finite(s) && s > 0) which((gaps - mean(tail)) / s > 6) else integer()
    k <- if (length(hits)) max(hits) else min(20L, length(f$d))
    k <- max(1L, min(k, length(f$d) - 1L))
  } else {
    if (length(rank) != 1L || !is.numeric(rank) || !is.finite(rank) || rank < 1)
      stop('rank must be auto or a positive integer', call. = FALSE)
    k <- min(as.integer(rank), kmax)
  }
  z <- sweep(f$u[, seq_len(min(20L, kmax)), drop = FALSE], 2L,
             f$d[seq_len(min(20L, kmax))], '*')
  rownames(z) <- colnames(x)
  list(U = sweep(f$u[, seq_len(k), drop = FALSE], 2L, f$d[seq_len(k)], '*'),
       V = f$v[, seq_len(k), drop = FALSE], embedding = z, rank = k)
}

# The full Walktrap merge history, including disconnected components, is kept.
# SuperCell 2.0 SCimplify_for_Seurat_v5.R uses floor(n / gamma) and cut_at.
.dk_lean_graph_membership <- function(graph, group, gamma = 150) {
  cells <- rownames(graph); n <- length(cells)
  if (is.null(cells) || !identical(cells, colnames(graph)) || anyDuplicated(cells))
    stop('graph must have identical unique row and column cell names', call. = FALSE)
  if (length(gamma) != 1L || !is.finite(gamma) || gamma < 2)
    stop('gamma must be >= 2', call. = FALSE)
  group <- .dk_align_vector(group, cells, 'group')
  if (is.null(group)) group <- rep('all', n)
  if (anyNA(group)) stop('broad groups cannot be missing', call. = FALSE)
  group <- as.character(group); names(group) <- cells
  mem <- integer(n); trees <- hierarchies <- graphs <- info <- list(); offset <- 0L
  for (s in unique(group)) {
    ids <- which(group == s)
    g <- igraph::graph_from_adjacency_matrix(graph[ids, ids, drop = FALSE],
                                            mode = 'undirected', weighted = TRUE, diag = FALSE)
    target <- max(1L, floor(length(ids) / gamma), igraph::components(g)$no)
    h <- if (igraph::ecount(g)) igraph::cluster_walktrap(g, steps = 4L, merges = TRUE) else NULL
    m <- if (is.null(h) || target >= length(ids)) seq_along(ids) else igraph::cut_at(h, no = target)
    m <- as.integer(factor(m, levels = unique(m)))
    mem[ids] <- offset + m; offset <- offset + max(m)
    hierarchies[s] <- list(h); trees[s] <- list(.dk_index_walktrap_tree(h, cells[ids]))
    graphs[s] <- list(g)
    info[[s]] <- data.frame(stratum=s, n_cells=length(ids), target_memberships=target,
                           observed_memberships=length(unique(m)), approximate=FALSE,
                           anchor_n=length(ids))
  }
  names(mem) <- cells
  tab <- as.data.frame(table(mem)); names(tab) <- c('membership','n_cells')
  tab$membership <- as.integer(as.character(tab$membership))
  structure(list(membership=mem, membership_table=tab, cell_stratum=group,
                 hierarchies=hierarchies, tree_indices=trees, graphs=graphs,
                 strata=do.call(rbind,info), settings=list(gamma=gamma, method='walktrap',
                   approximate=FALSE, has_hard_stratum=length(unique(group))>1L)),
            class='DropoutKillerMembership')
}

.dk_lean_rna_geometry <- function(z, group, gamma = 150, k_knn = 5L) {
  if (is.null(group)) group <- rep('all', nrow(z))
  if (anyNA(group)) stop('broad groups cannot be missing',call.=FALSE)
  edge <- list()
  for (s in unique(group)) {
    ids <- which(group == s)
    g <- .dk_knn_graph(z[ids,,drop=FALSE], k_knn)
    e <- igraph::as_edgelist(g, names=FALSE)
    if (nrow(e)) edge[[length(edge)+1L]] <- cbind(ids[e[,1]],ids[e[,2]])
  }
  e <- do.call(rbind,edge)
  if (is.null(e)) e <- matrix(integer(),0,2)
  a <- Matrix::sparseMatrix(i=c(e[,1],e[,2]),j=c(e[,2],e[,1]),x=1,
                           dims=c(nrow(z),nrow(z)),dimnames=list(rownames(z),rownames(z)))
  a@x[] <- 1
  list(embedding=z, membership_fit=.dk_lean_graph_membership(a,group,gamma), affinity=NULL)
}

# Own metacell is exhausted first; only then walk up to the nearest metacells.
# RNA uses Euclidean distance inside each eligible membership. WNN uses graph
# affinity and never replaces missing WNN links with RNA-only distances.
.dk_lean_neighbors <- function(geometry, k = 30L) {
  z <- geometry$embedding; fit <- geometry$membership_fit
  cells <- rownames(z); mem <- fit$membership[cells]; grp <- fit$cell_stratum[cells]
  if (anyNA(mem) || anyNA(grp)) stop('geometry cell alignment failed', call. = FALSE)
  if (length(k)!=1L || !is.finite(k) || k < 1) stop('neighbor k must be >= 1',call.=FALSE)
  affinity <- geometry$affinity; edges <- list(); borrowed <- integer(length(cells))
  for (s in unique(grp)) {
    ids <- which(grp==s); mids <- unique(mem[ids]); tree <- fit$tree_indices[[s]]
    pools <- lapply(mids,function(m)ids[mem[ids]==m])
    reps <- vapply(pools, function(p) p[1L], integer(1))
    centers <- t(vapply(pools,function(p)colMeans(z[p,,drop=FALSE]),numeric(ncol(z))))
    for (a in seq_along(mids)) {
      own <- pools[[a]]; others <- integer()
      if (length(own)-1L<k && length(mids)>1L && !is.null(tree)) {
        td <- .dk_tree_distance(tree,rep(cells[reps[a]],length(reps)),cells[reps])
        # Disconnected roots have no shared ancestor: do not invent proximity.
        connected <- vapply(reps,function(r)length(intersect(tree$ancestors[[cells[r]]],
                                     tree$ancestors[[cells[reps[a]]]]))>0L,logical(1))
        td[!connected] <- Inf; td[a] <- Inf
        tie <- if (is.null(affinity)) rowSums(sweep(centers,2,centers[a,],'-')^2) else
          -vapply(pools,function(p)sum(affinity[own,p,drop=FALSE]),numeric(1))
        others <- order(td,tie,mids); others <- others[is.finite(td[others])]
      }
      for (c in own) {
        rank_pool <- function(pool) {
          pool <- setdiff(pool,c)
          if (!length(pool)) return(pool)
          if (is.null(affinity)) {
            d <- rowSums(sweep(z[pool,,drop=FALSE],2,z[c,],'-')^2)
            pool[order(d,pool)]
          } else {
            w <- as.numeric(affinity[c,pool]); valid <- which(w>0)
            pool[valid[order(-w[valid],pool[valid])]]
          }
        }
        take <- head(rank_pool(own),k); same <- length(take)
        # In WNN an own membership can contain many cells but few usable links.
        # Such sparse support must also allow hierarchy borrowing.
        if (!is.null(affinity) && length(take)<k && !length(others) && length(mids)>1L && !is.null(tree)) {
          td <- .dk_tree_distance(tree,rep(cells[c],length(reps)),cells[reps]);td[a]<-Inf
          connected <- vapply(reps,function(r)length(intersect(tree$ancestors[[cells[c]]],tree$ancestors[[cells[r]]]))>0L,logical(1))
          td[!connected]<-Inf
          others <- order(td,mids);others<-others[is.finite(td[others])]
        }
        for (b in others) {
          if (length(take)>=k) break
          take <- c(take,head(rank_pool(pools[[b]]),k-length(take)))
        }
        if (!length(take)) next
        if (is.null(affinity)) {
          ds <- sqrt(rowSums(sweep(z[take,,drop=FALSE],2,z[c,],'-')^2));h<-stats::median(ds)
          w <- if(h>0)exp(-ds^2/(2*h*h))else as.numeric(ds==0)
        } else w <- as.numeric(affinity[c,take])
        edges[[length(edges)+1L]] <- cbind(i=c,j=take,w=w/sum(w))
        borrowed[c] <- length(take)-same
      }
    }
  }
  e<-do.call(rbind,edges)
  if(is.null(e))e<-matrix(numeric(),0,3)
  list(W=Matrix::sparseMatrix(i=e[,1],j=e[,2],x=e[,3],dims=c(length(cells),length(cells))),
       borrowed=borrowed)
}

.dk_lean_context <- function(x, geometry) {
  mem <- as.integer(factor(geometry$membership_fit$membership[colnames(x)]))
  grp <- as.integer(factor(geometry$membership_fit$cell_stratum[colnames(x)]))
  mn<-tabulate(mem);gn<-tabulate(grp);n<-ncol(x)
  pm<-(x>0)%*%Matrix::sparseMatrix(i=seq_len(n),j=mem,x=1)
  pg<-(x>0)%*%Matrix::sparseMatrix(i=seq_len(n),j=grp,x=1)
  list(mem=mem,grp=grp,mn=mn,gn=gn,pm=pm,pg=pg)
}

.dk_lean_features <- function(x, coordinates, geometry, svd_fit, neighbors=NULL, context=NULL) {
  if (is.null(neighbors)) neighbors <- .dk_lean_neighbors(geometry)
  if (is.null(context)) context <- .dk_lean_context(x,geometry)
  mem<-context$mem;grp<-context$grp;mn<-context$mn;gn<-context$gn;pm<-context$pm;pg<-context$pg
  out<-data.frame(z=numeric(nrow(coordinates)),local_hierarchy=numeric(nrow(coordinates)),H=numeric(nrow(coordinates)))
  if(!nrow(coordinates))return(out)
  ord<-order(coordinates[,1],method='radix');sorted<-coordinates[ord,1]
  ends<-c(which(diff(sorted)!=0),length(sorted));starts<-c(1L,head(ends,-1L)+1L)
  missing_neighbors<-Matrix::rowSums(neighbors$W)==0
  for (b in seq_along(starts)) {
    ix<-ord[seq.int(starts[b],ends[b])];g<-coordinates[ix[1L],1];cc<-coordinates[ix,2]
    lr<-as.numeric(svd_fit$U%*%svd_fit$V[g,])
    tau<-abs(stats::quantile(lr,.001,names=FALSE));sig<-stats::sd(lr[lr<0])
    if(!is.finite(sig)||sig<=0)sig<-stats::sd(lr)
    if(!is.finite(sig)||sig<1e-8)sig<-1e-8
    out$z[ix]<-(lr[cc]-tau)/sig
    support<-as.numeric(neighbors$W%*%as.numeric(x[g,]>0))
    support[missing_neighbors]<-.5
    out$local_hierarchy[ix]<-stats::qlogis((pmin(1,pmax(0,support[cc]))+.01)/1.02)
    a<-(as.numeric(pm[g,mem[cc]])+.5)/(mn[mem[cc]]+1)
    b<-(as.numeric(pg[g,grp[cc]])+.5)/(gn[grp[cc]]+1)
    out$H[ix]<-mn[mem[cc]]/(mn[mem[cc]]+50)*(stats::qlogis(a)-stats::qlogis(b))
  }
  out
}

#' Fit the balanced Lean hierarchy detector and a proxy-FPR threshold
#' @param features Data frame with z, local_hierarchy and H benchmark features.
#' @param truth Logical labels: artificially lost positives versus proxy negatives.
#' @param fpr Training proxy false-positive target, strictly between zero and one.
#' @param method RNA or WNN detector name; models are not interchangeable.
#' @return A reusable DropoutKillerLeanModel. Its score is not a biological posterior.
#' @export
fit_lean_detector <- function(features, truth, fpr=.01, method='Supercell_hierarchy_Lean_membership') {
  cols<-c('z','local_hierarchy','H');X<-as.matrix(features[,cols,drop=FALSE])
  if(!is.logical(truth)||anyNA(truth)||length(truth)!=nrow(X)||!any(truth)||all(truth)||any(!is.finite(X)))
    stop('finite features and both logical truth classes are required',call.=FALSE)
  if(length(fpr)!=1L||!is.finite(fpr)||fpr<=0||fpr>=1)stop('fpr must be in (0,1)',call.=FALSE)
  if(!method%in%c(.dk_lean_method(),.dk_lean_method(TRUE)))stop('invalid Lean method',call.=FALSE)
  w<-ifelse(truth,.5/sum(truth),.5/sum(!truth));mu<-colSums(X*w)
  scale<-sqrt(pmax(colSums(sweep(X,2,mu,'-')^2*w),1e-16))
  # Constant features contribute zero; dropping them prevents a singular fit.
  active<-scale>1e-8
  design<-cbind(1,sweep(sweep(X[,active,drop=FALSE],2,mu[active],'-'),2,scale[active],'/'))
  fit<-stats::glm.fit(design,as.numeric(truth),weights=w*length(truth),
                     family=stats::quasibinomial(),control=stats::glm.control(maxit=100,epsilon=1e-9))
  if(!isTRUE(fit$converged)||any(!is.finite(fit$coefficients)))stop('Lean calibration did not converge',call.=FALSE)
  beta<-stats::setNames(numeric(4),c('(Intercept)',cols));beta[-1][active]<-fit$coefficients[-1]/scale[active]
  beta[1]<-fit$coefficients[1]-sum(beta[-1]*mu)
  scores<-as.numeric(beta[1]+X%*%beta[-1]);neg<-sort(scores[!truth])
  cutoff<-neg[max(1L,ceiling((1-fpr)*length(neg)))]
  structure(list(beta=beta,threshold=cutoff,fpr=fpr,method=method,converged=TRUE,
                 train_counts=c(positive=sum(truth),negative=sum(!truth)),
                 feature_version='lean-hierarchy-v1',scope='training proxy-FPR, not biological FDR'),
            class='DropoutKillerLeanModel')
}

.dk_validate_lean_model <- function(model, method) {
  if(!inherits(model,'DropoutKillerLeanModel')||!isTRUE(model$converged)||
     !identical(model$method,method)||!identical(model$feature_version,'lean-hierarchy-v1')||
     !identical(names(model$beta),c('(Intercept)','z','local_hierarchy','H'))||
     any(!is.finite(model$beta))||length(model$threshold)!=1L||!is.finite(model$threshold))
    stop('supply a converged fit_lean_detector model for this RNA/WNN method',call.=FALSE)
  model
}

.dk_validate_lean_geometry <- function(geometry,cells,group=NULL) {
  fit<-geometry$membership_fit
  if(!inherits(fit,'DropoutKillerMembership')||
     !identical(rownames(geometry$embedding),cells)||
     !identical(names(fit$membership),cells)||!identical(names(fit$cell_stratum),cells)||
     anyNA(fit$membership)||anyNA(fit$cell_stratum))stop('geometry cell order or membership alignment failed',call.=FALSE)
  if(!is.null(group)) {
    group<-.dk_align_vector(group,cells,'group')
    if(!identical(as.character(group),unname(as.character(fit$cell_stratum))))
      stop('geometry broad groups differ from supplied group',call.=FALSE)
  }
  if(any(vapply(split(fit$cell_stratum,fit$membership),function(s)length(unique(s))!=1L,logical(1))))
    stop('membership crosses broad groups',call.=FALSE)
  a<-geometry$affinity
  if(!is.null(a)) {
    if(!identical(rownames(a),cells)||!identical(colnames(a),cells)||!Matrix::isSymmetric(a))
      stop('WNN affinity cell order or symmetry invalid',call.=FALSE)
    e<-Matrix::summary(a)
    if(any(!is.finite(e$x)|e$x<0)||any(e$i==e$j & e$x!=0)||
       any(fit$cell_stratum[e$i]!=fit$cell_stratum[e$j] & e$x!=0))
      stop('WNN affinity contains invalid, self, or cross-group edges',call.=FALSE)
  }
  invisible(geometry)
}

#' Detect recoverable zeros with SuperCell hierarchy and Lean membership features
#' @param x Nonnegative gene-by-cell expression; raw counts if normalize is TRUE.
#' @param embedding RNA cell embedding, or NULL to derive it from the ALRA SVD.
#' @param group Broad cell classes. Fine subclusters are not inferred as boundaries.
#' @param model A model returned by fit_lean_detector.
#' @param geometry Optional precomputed internal geometry, including WNN affinity.
#' @param normalize Apply library-size normalization and log1p.
#' @param gamma Cells per metacell, with floor(n/gamma) Walktrap cut.
#' @param k_knn Graph neighbor count for RNA memberships.
#' @param neighbor_k Maximum support neighbors, borrowing only on shortage.
#' @param rank ALRA rank or auto.
#' @param seed Random seed for truncated SVD.
#' @return DropoutKillerDetection with called zero events and fitted model provenance.
#' @export
supercell_lean_detect <- function(x, embedding=NULL, group=NULL, model=NULL,
                                  geometry=NULL, normalize=TRUE, gamma=150,
                                  k_knn=5L, neighbor_k=30L, rank='auto', seed=1L) {
  x<-.dk_validate_expression(x);nm<-.dk_names(x);dimnames(x)<-list(nm$genes,nm$cells)
  group<-.dk_align_vector(group,nm$cells,'group')
  if(normalize)x<-.dk_alra_library_log(x)
  sv<-.dk_lean_svd(x,rank,seed)
  if(is.null(geometry)) {
    z<-if(is.null(embedding))sv$embedding else .dk_align_embedding(embedding,nm$cells)
    geometry<-.dk_lean_rna_geometry(z,group,gamma,k_knn)
  }
  .dk_validate_lean_geometry(geometry,nm$cells,group)
  method<-.dk_lean_method(!is.null(geometry$affinity));model<-.dk_validate_lean_model(model,method)
  nb<-.dk_lean_neighbors(geometry,neighbor_k);events<-list();context<-.dk_lean_context(x,geometry)
  # Bound the reconstruction/event working set by one gene; never densify all zeros.
  for(g in seq_len(nrow(x))) {
    cc<-which(as.numeric(x[g,])==0);if(!length(cc))next
    ij<-cbind(g,cc);f<-.dk_lean_features(x,ij,geometry,sv,nb,context)
    s<-as.numeric(model$beta[1]+as.matrix(f)%*%model$beta[-1]);keep<-which(s>model$threshold)
    if(!length(keep))next
    events[[length(events)+1L]]<-data.frame(i=g,j=cc[keep],gene=rownames(x)[g],cell=colnames(x)[cc[keep]],
      membership=unname(geometry$membership_fit$membership[cc[keep]]),
      lowrank=as.numeric(sv$U[cc[keep],,drop=FALSE]%*%sv$V[g,]),threshold=model$threshold,
      null_sigma=NA_real_,z_score=f$z[keep],p_value=NA_real_,q_value=NA_real_,
      confidence=stats::plogis(s[keep]),confidence_fallback=FALSE,variance_weight=NA_real_,
      lean_score=s[keep],local_hierarchy=f$local_hierarchy[keep],H=f$H[keep])
  }
  ev<-do.call(rbind,events)
  if(is.null(ev))ev<-data.frame(i=integer(),j=integer(),gene=character(),cell=character(),membership=integer(),
    lowrank=numeric(),threshold=numeric(),null_sigma=numeric(),z_score=numeric(),p_value=numeric(),q_value=numeric(),
    confidence=numeric(),confidence_fallback=logical(),variance_weight=numeric(),lean_score=numeric(),local_hierarchy=numeric(),H=numeric())
  structure(list(events=ev,dimensions=dim(x),dimnames=dimnames(x),model=model,
    membership_fit=geometry$membership_fit,membership_stats=geometry$membership_fit$strata,
    settings=list(detection_method=method,detection_scope='hierarchy_within_broad_group',
                  rank=sv$rank,neighbor_k=neighbor_k,borrowed_cells=sum(nb$borrowed>0),
                  score_type='balanced_logistic_score_not_posterior',fpr=model$fpr)),
    class='DropoutKillerDetection')
}
