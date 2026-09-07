#!/usr/bin/env python3
import argparse, os, json, time, glob, warnings
import numpy as np
import pandas as pd
import h5py
from scipy import sparse
from scipy.special import digamma, expit, logit
from scipy.optimize import root
from scipy.stats import gamma, norm
from sklearn.utils.extmath import randomized_svd
from sklearn.neighbors import NearestNeighbors
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
warnings.filterwarnings('ignore')
POINT=float(np.log(1.01))
RHO_LIST=[0.85,0.70,0.55]
LABEL_ORDER=['B_cell','Dendritic_cell','HSPC','Monocyte','NK_cell','T_cell']
PASS_RULE={'tpr_fpr1_min_delta':0.02,'tpr_fpr5_min_delta':0.01,'auprc_residual_reduction_min':0.20,'auroc_residual_reduction_min':0.20,'paired_wins_min':24,'default_fpr_margin_vs_alra':0.005,'require_default_tpr_ge_alra':True,'require_each_rho_positive':True}
def lg(x): return logit(np.clip(x,1e-6,1-1e-6))
def dec(a): return np.asarray([x.decode() if isinstance(x,(bytes,np.bytes_)) else str(x) for x in a])
def read_input(inp):
    meta=pd.read_csv(os.path.join(inp,'metadata.csv'))
    emb=pd.read_csv(os.path.join(inp,'embedding20.csv'),index_col=0).loc[meta.cell.values].values.astype(np.float32)
    mem=pd.read_csv(os.path.join(inp,'membership.csv')).set_index('cell').loc[meta.cell.values,'membership'].values.astype(np.int32)
    with h5py.File(os.path.join(inp,'pbmc_full.h5'),'r') as h:
        g=h['matrix']; shape=tuple(g['shape'][...]); data=g['data'][...]; indices=g['indices'][...]; indptr=g['indptr'][...]; ft=dec(g['features']['feature_type'][...]); names=dec(g['features']['name'][...]); bar=dec(g['barcodes'][...])
    M=sparse.csc_matrix((data,indices,indptr),shape=shape); gi=np.where(ft=='Gene Expression')[0]; G=M[gi,:]; genes=names[gi]; bmap={b:i for i,b in enumerate(bar)}
    missing=[b for b in meta.cell if b not in bmap]
    if missing: raise RuntimeError(f'{len(missing)} metadata barcodes absent from 10x matrix; first={missing[:3]}')
    ci=np.array([bmap[b] for b in meta.cell],dtype=np.int32); C=G[:,ci].tocsc().astype(np.int32); lab=meta.cell_type.values.astype(str); lum={x:i for i,x in enumerate(LABEL_ORDER)}
    if set(lab)!=set(LABEL_ORDER): raise RuntimeError(f'label mismatch {sorted(set(lab))}')
    cli=np.array([lum[x] for x in lab],dtype=np.int16); return C,genes,meta,lab,cli,emb,mem
def binary_matrix(C):
    B=C.copy().astype(np.int8); B.data[:]=1; return B
def provenance(C,lab):
    G,N=C.shape; B=binary_matrix(C); det=np.asarray(B.sum(1)).ravel()/N; keep=det>=0.05; sub=B[keep]; zero=1-sub.nnz/(sub.shape[0]*sub.shape[1]); counts={x:int(np.sum(lab==x)) for x in LABEL_ORDER}
    checks={'n_genes_all':G,'n_cells':N,'n_genes_detection_ge_5pct':int(keep.sum()),'zero_fraction_detection_ge_5pct':float(zero),'cell_type_counts':counts,'expected_dim_pass':bool(G==36601 and N==10412),'expected_old_filter_n_pass':bool(int(keep.sum())==8793),'expected_old_zero_pass':bool(abs(zero-0.7879)<5e-4),'expected_cell_counts_pass':bool(counts=={'B_cell':884,'Dendritic_cell':304,'HSPC':26,'Monocyte':3326,'NK_cell':468,'T_cell':5404})}
    checks['all_pass']=all(checks[k] for k in ['expected_dim_pass','expected_old_filter_n_pass','expected_old_zero_pass','expected_cell_counts_pass']); return checks,det
def supports(C,cli,L=6):
    B=binary_matrix(C); N=C.shape[1]; one=sparse.csr_matrix((np.ones(N,np.int8),(np.arange(N),cli)),shape=(N,L)); return np.asarray((B@one).toarray(),dtype=np.int32)
def build_neg_controls(C0,cli,npos0,cap=100,seed=20260907):
    G,N=C0.shape; L=npos0.shape[1]; nlab=np.bincount(cli,minlength=L); D=npos0/nlab[None,:]; Cr=C0.tocsr(); gs=[];cs=[];ls=[];rng=np.random.default_rng(seed)
    for li in range(L):
        other=np.max(np.delete(D,li,axis=1),axis=1); qual=np.where((D[:,li]<=.005)&(other>=.20))[0]; cells=np.where(cli==li)[0]
        for g in qual:
            pos=Cr.indices[Cr.indptr[g]:Cr.indptr[g+1]]; z=cells[~np.isin(cells,pos,assume_unique=False)]
            if len(z)>cap: z=rng.choice(z,cap,replace=False)
            if len(z): gs.append(np.full(len(z),g,np.int32));cs.append(z.astype(np.int32));ls.append(np.full(len(z),li,np.int8))
    return (np.concatenate(gs),np.concatenate(cs),np.concatenate(ls)) if gs else (np.array([],np.int32),np.array([],np.int32),np.array([],np.int8))
def make_mask(C0,rho,seed,cli,npos0):
    rng=np.random.default_rng(seed);G,N=C0.shape;L=npos0.shape[1];coo=C0.tocoo(copy=False);g=coo.row.astype(np.int32);c=coo.col.astype(np.int32);v=coo.data.astype(np.int32);elig=npos0[g,cli[c]]>=5;p=np.zeros(len(v),np.float32);p[elig]=np.power(1-rho,v[elig].astype(np.float32));sel=(rng.random(len(v))<p)&elig;poscell=np.diff(C0.indptr);ix=np.flatnonzero(sel)
    if len(ix):
        order=np.argsort(c[ix],kind='mergesort');sx=ix[order];sc=c[sx];starts=np.r_[0,np.flatnonzero(np.diff(sc))+1];ends=np.r_[starts[1:],len(sx)]
        for a,b in zip(starts,ends):
            cell=int(sc[a]);ids=sx[a:b];cap=int(np.floor(.20*poscell[cell]))
            if len(ids)>cap: sel[ids]=False; sel[rng.choice(ids,cap,replace=False) if cap>0 else np.array([],int)]=True
    ix=np.flatnonzero(sel)
    if len(ix):
        keys=g[ix].astype(np.int64)*L+cli[c[ix]];order=np.argsort(keys,kind='mergesort');sx=ix[order];sk=keys[order];starts=np.r_[0,np.flatnonzero(np.diff(sk))+1];ends=np.r_[starts[1:],len(sx)]
        for a,b in zip(starts,ends):
            key=int(sk[a]);gg=key//L;li=key%L;ids=sx[a:b];maxmask=max(0,int(npos0[gg,li])-3)
            if len(ids)>maxmask: sel[ids]=False; sel[rng.choice(ids,maxmask,replace=False) if maxmask>0 else np.array([],int)]=True
    capall=int(np.floor(.05*G*N));ix=np.flatnonzero(sel)
    if len(ix)>capall: sel[ix]=False;sel[rng.choice(ix,capall,replace=False)]=True
    mg=g[sel];mc=c[sel];mv=v[sel];keep=~sel;C=sparse.csc_matrix((v[keep],(g[keep],c[keep])),shape=C0.shape,dtype=np.int32);return C,mg,mc,mv
def sc_fit_sparse(posvals,n):
    npos=len(posvals);nzero=n-npos;xp=np.log(1.01+posvals.astype(np.float64));meanx=(nzero*POINT+xp.sum())/n
    if abs(meanx-POINT)<1e-2:return np.nan,np.nan,True
    rate=nzero/n
    if rate>0.95 or npos==0:return np.nan,np.nan,True
    if rate==0:rate=.01
    alpha,beta=1.5,1.;mu=float(xp.mean());std=float(xp.std()) or .01;x=np.r_[POINT,xp];mult=np.r_[float(nzero),np.ones(npos)];eps=10.;old=0.;it=0
    while eps>.5:
        p1=rate*gamma.pdf(x,a=alpha,scale=1/beta);p2=(1-rate)*norm.pdf(x,mu,std);den=p1+p2;z=np.divide(p1,den,out=np.zeros_like(p1),where=(den>0)&np.isfinite(den));z[p1==0]=0;w0=mult*z;w1=mult*(1-z);s0=w0.sum();s1=w1.sum()
        if s0<=0 or s1<=0:return np.nan,np.nan,True
        rate2=s0/n;mu2=float(np.sum(w1*x)/s1);std2=float(np.sqrt(np.sum(w1*(x-mu2)**2)/s1));tt=float(np.sum(w0*x));tu=float(np.sum(w0*np.log(x)))
        if tt<=0 or not np.isfinite(std2):return np.nan,np.nan,True
        tv=-tu/s0-np.log(s0/tt)
        if tv<=0:a2=20.
        else:
            a0=(3-tv+np.sqrt((tv-3)**2+24*tv))/(12*tv)
            if a0>=20:a2=20.
            else:a2=float(root(lambda a:np.log(a)-digamma(a)-tv,np.array([.9,1.1])*a0).x[0])
        b2=s0/tt*a2;rate,alpha,beta,mu,std=rate2,a2,b2,mu2,max(std2,1e-8);new=rate*gamma.pdf(x,a=alpha,scale=1/beta)+(1-rate)*norm.pdf(x,mu,std)
        if np.any(new<=0)|np.any(~np.isfinite(new)):return np.nan,np.nan,True
        ll=float(np.sum(mult*np.log10(new)));eps=(ll-old)**2;old=ll;it+=1
        if it>100:break
    ga=rate*gamma.pdf(POINT,a=alpha,scale=1/beta);no=(1-rate)*norm.pdf(POINT,mu,std);d=ga/(ga+no) if ga+no>0 else np.nan;return float(d),float(rate),False
def scgacl_params(C,cli):
    G,N=C.shape;L=6;d=np.full((G,L),np.nan,np.float32);lam=np.full((G,L),np.nan,np.float32);invalid=np.ones((G,L),bool)
    for li in range(L):
        ids=np.where(cli==li)[0];R=C[:,ids].tocsr();n=len(ids)
        for g in range(G):
            vals=R.data[R.indptr[g]:R.indptr[g+1]];dd,ll,ii=sc_fit_sparse(vals,n);d[g,li]=dd;lam[g,li]=ll;invalid[g,li]=ii
    R=C.tocsr();dg=np.full(G,np.nan,np.float32);glam=np.full(G,np.nan,np.float32);invg=np.ones(G,bool)
    for g in range(G):
        vals=R.data[R.indptr[g]:R.indptr[g+1]];dd,ll,ii=sc_fit_sparse(vals,N);dg[g]=dd;glam[g]=ll;invg[g]=ii
    return d,lam,dg,glam,invalid,invg
def normalize_sparse(C,lib):
    X=C.T.tocsr().astype(np.float32);rows=np.repeat(np.arange(X.shape[0]),np.diff(X.indptr));X.data=np.log1p(X.data*10000.0/np.maximum(lib[rows],1.0)).astype(np.float32);return X
def choose_k(X,seed):
    K=min(100,min(X.shape)-1);_,S,_=randomized_svd(X,n_components=K,n_iter=2,random_state=seed);dif=S[:-1]-S[1:];ns=min(80,K-4);idx=np.arange(max(1,ns-1),len(dif));mu=float(dif[idx].mean());sd=float(dif[idx].std(ddof=1)) if len(idx)>1 else 0
    if not np.isfinite(sd) or sd<=0:return 1
    z=(dif-mu)/sd;hit=np.where(z>6)[0];return int(hit.max()+1 if len(hit) else 1)
def alra_events(C,lib,eg,ec,seed,block=384):
    X=normalize_sparse(C,lib);k=choose_k(X,seed);U,S,Vt=randomized_svd(X,n_components=k,n_iter=10,random_state=seed+100000);US=(U*S).astype(np.float32);Vt=Vt.astype(np.float32);Xc=X.tocsc();G=X.shape[1];outz=np.empty(len(eg),np.float32);outcall=np.zeros(len(eg),bool);order=np.argsort(eg,kind='mergesort');sg=eg[order]
    for a in range(0,G,block):
        b=min(G,a+block);lo=np.searchsorted(sg,a,'left');hi=np.searchsorted(sg,b,'left')
        if hi<=lo:continue
        ei=order[lo:hi];loc=eg[ei]-a;LR=US@Vt[:,a:b];tau=np.abs(np.quantile(LR,.001,axis=0));neg=np.where(LR<0,LR,np.nan);negsd=np.nanstd(neg,axis=0,ddof=1);fb=np.nanmedian(negsd[np.isfinite(negsd)&(negsd>0)]);negsd=np.where(np.isfinite(negsd)&(negsd>0),negsd,fb if np.isfinite(fb) else 1.);lev=LR[ec[ei],loc];outz[ei]=((lev-tau[loc])/negsd[loc]).astype(np.float32);pos=LR>tau[None,:];cnt=pos.sum(0);s1=np.where(pos,LR,0).sum(0);s2=np.where(pos,LR*LR,0).sum(0);mu1=np.divide(s1,cnt,out=np.full(b-a,np.nan),where=cnt>0);var1=np.divide(s2-s1*s1/np.maximum(cnt,1),cnt-1,out=np.full(b-a,np.nan),where=cnt>1);sd1=np.sqrt(np.maximum(var1,0));mu2=np.full(b-a,np.nan);sd2=np.full(b-a,np.nan)
        for j in range(a,b):
            vals=Xc.data[Xc.indptr[j]:Xc.indptr[j+1]].astype(float);q=j-a
            if len(vals):mu2[q]=vals.mean()
            if len(vals)>1:sd2[q]=vals.std(ddof=1)
        good=np.isfinite(sd1)&np.isfinite(sd2)&(sd1!=0);ratio=np.ones(b-a);add=np.zeros(b-a);ratio[good]=sd2[good]/sd1[good];add[good]=-mu1[good]*ratio[good]+mu2[good];scaled=lev*ratio[loc]+add[loc];outcall[ei]=(lev>tau[loc])&(scaled>0)
    return outz,outcall,k
def build_W(emb,cli,k=30):
    N=len(cli);rows=[];cols=[];dat=[]
    for li in range(6):
        ids=np.where(cli==li)[0];kk=min(k+1,len(ids));nn=NearestNeighbors(n_neighbors=kk).fit(emb[ids]);dist,ind=nn.kneighbors(emb[ids])
        for aa,c in enumerate(ids):
            js=ids[ind[aa]];ds=dist[aa];keep=js!=c;js=js[keep];ds=ds[keep]
            if not len(js):continue
            h=np.median(ds);h=h if h>0 else max(np.mean(ds),1e-6);w=np.exp(-(ds*ds)/(2*h*h));w/=w.sum();rows.extend([c]*len(js));cols.extend(js.tolist());dat.extend(w.tolist())
    return sparse.csr_matrix((dat,(rows,cols)),shape=(N,N),dtype=np.float32)
def neighbor_events(C,W,eg,ec,block=384):
    B=binary_matrix(C).T.tocsr();out=np.empty(len(eg),np.float32);order=np.argsort(eg,kind='mergesort');sg=eg[order];G=C.shape[0]
    for a in range(0,G,block):
        b=min(G,a+block);lo=np.searchsorted(sg,a);hi=np.searchsorted(sg,b)
        if hi<=lo:continue
        ei=order[lo:hi];S=(W@B[:,a:b]).toarray();val=S[ec[ei],eg[ei]-a];out[ei]=lg((val+.01)/1.02).astype(np.float32)
    return out
def context_events(C,lib,cli,mem,scd,scl,dg,glam,eg,ec):
    G,N=C.shape;B=binary_matrix(C);L=6;nlab=np.bincount(cli,minlength=L);labone=sparse.csr_matrix((np.ones(N,np.int8),(np.arange(N),cli)),shape=(N,L));npl=np.asarray((B@labone).toarray(),dtype=np.int32);alln=np.asarray(B.sum(1)).ravel();M=int(mem.max());mz=mem-1;msize=np.bincount(mz,minlength=M);memone=sparse.csr_matrix((np.ones(N,np.int8),(np.arange(N),mz)),shape=(N,M));npm=np.asarray((B@memone).toarray(),dtype=np.int32);dsh=np.empty((G,L),np.float32);ash=np.empty((G,L),np.float32)
    for li in range(L):
        r=npl[:,li]/(npl[:,li]+20.);dl=np.where(np.isfinite(scd[:,li]),scd[:,li],dg);dl=np.where(np.isfinite(dl),dl,.01);gd=np.where(np.isfinite(dg),dg,.01);lm=np.where(np.isfinite(scl[:,li]),scl[:,li],glam);lm=np.where(np.isfinite(lm),lm,.5);gg=np.where(np.isfinite(glam),glam,.5);dsh[:,li]=expit(r*lg(dl)+(1-r)*lg(gd));ash[:,li]=expit(r*lg(1-lm)+(1-r)*lg(1-gg))
    el=cli[ec];em=mz[ec];rr=msize[em]/(msize[em]+50.);leaf=(npm[eg,em]+.5)/(msize[em]+1.);pp=(npl[eg,el]+.5)/(nlab[el]+1.);delta=rr*(lg(leaf)-lg(pp));hd=lg(dsh[eg,el])+delta;ha=lg(ash[eg,el])+delta;theta=np.zeros((G,L),np.float64);alpha=np.zeros((G,L),np.float64)
    for li in range(L):
        ids=np.where(cli==li)[0];R=C[:,ids].tocsr().astype(np.float64);Lv=lib[ids].astype(float);sl=Lv.sum();sl2=np.dot(Lv,Lv);sy=np.asarray(R.sum(1)).ravel();R2=R.copy();R2.data**=2;sy2=np.asarray(R2.sum(1)).ravel();syL=np.asarray(R@Lv).ravel();th=sy/max(sl,1.);den=th*th*sl2;num=sy2-sy-2*th*syL+den;al=np.maximum(np.divide(num,den+1e-12),0);theta[:,li]=th;alpha[:,li]=al
    mu=theta[eg,el]*lib[ec];al=alpha[eg,el];q=np.empty(len(eg),float);po=al<1e-8;q[po]=np.exp(-mu[po]);aa=al[~po];q[~po]=np.power(1+aa*mu[~po],-1/aa);return hd.astype(np.float32),ha.astype(np.float32),np.log(np.clip(q,1e-12,1)).astype(np.float32)
def event_features(inp,rho,seed,neg_cap=100):
    t=time.time();C0,genes,meta,lab,cli,emb,mem=read_input(inp);prov,det=provenance(C0,lab)
    if not prov['all_pass']:raise RuntimeError('PROVENANCE FAILURE: '+json.dumps(prov))
    npos0=supports(C0,cli);ng,nc,nl=build_neg_controls(C0,cli,npos0,cap=neg_cap);C,mg,mc,mv=make_mask(C0,rho,seed,cli,npos0);lost=np.bincount(mc,weights=mv,minlength=C0.shape[1]);lib=np.maximum(np.asarray(C0.sum(0)).ravel()-lost,1);eg=np.r_[mg,ng].astype(np.int32);ec=np.r_[mc,nc].astype(np.int32);y=np.r_[np.ones(len(mg),np.int8),np.zeros(len(ng),np.int8)];orig=np.r_[mv.astype(np.int32),np.zeros(len(ng),np.int32)];scd,scl,dg,glam,inv,invg=scgacl_params(C,cli);az,acall,rank=alra_events(C,lib,eg,ec,seed);W=build_W(emb,cli,30);neigh=neighbor_events(C,W,eg,ec);hd,ha,lq=context_events(C,lib,cli,mem,scd,scl,dg,glam,eg,ec);el=cli[ec];rawsc=scd[eg,el];scscore=np.where(np.isfinite(rawsc),rawsc,0).astype(np.float32);sccall=scscore>.5;diag={'rho':rho,'seed':seed,'n_genes_all':C0.shape[0],'n_cells':C0.shape[1],'masked_n':len(mg),'negative_n':len(ng),'alra_rank_auto':rank,'baseline_zero_fraction_all':float(1-C0.nnz/(C0.shape[0]*C0.shape[1])),'masked_zero_fraction_all':float(1-C.nnz/(C.shape[0]*C.shape[1])),'scgacl_invalid_gene_lineage_fraction':float(inv.mean()),'scgacl_invalid_global_gene_fraction':float(invg.mean()),'runtime_sec':time.time()-t,'provenance':prov,'neg_cap_per_gene_lineage':neg_cap};return {'y':y,'g':eg,'c':ec,'orig_count':orig,'alra_z':az,'alra_call':acall,'sc_score':scscore,'sc_call':sccall,'hier_d':hd,'hier_a':ha,'neigh':neigh,'logq':lq,'det':det,'mem':mem,'cli':cli,'diag':diag}
def sample_cal(f,seed,n=10000):
    y=f['y'];rng=np.random.default_rng(seed);pp=np.where(y==1)[0];nn=np.where(y==0)[0];k=min(n,len(pp),len(nn));ix=np.r_[rng.choice(pp,k,False),rng.choice(nn,k,False)];return y[ix],f['alra_z'][ix],f['neigh'][ix]
def calc_metrics(y,score,call):
    p=score[y==1];n=score[y==0];o={'auroc':roc_auc_score(y,score),'auprc':average_precision_score(y,score),'tpr_default':float(call[y==1].mean()),'fpr_default':float(call[y==0].mean()),'precision_default':float(np.sum(call&(y==1))/max(np.sum(call),1))}
    for q in [.01,.05]:th=np.quantile(n,1-q,method='higher');o[f'th_fpr{int(q*100)}']=float(th);o[f'tpr_fpr{int(q*100)}']=float(np.mean(p>th))
    return o
def save_feature_mode(args):
    f=event_features(args.input,args.rho,args.seed,args.neg_cap);y,a,n=sample_cal(f,args.seed+20260907,args.cal_per_class);np.savez_compressed(args.output,y=y,alra_z=a,neigh=n,rho=args.rho,seed=args.seed);json.dump(f['diag'],open(args.output+'.diag.json','w'),indent=2)
def calibrate(args):
    files=glob.glob(os.path.join(args.cal_dir,'**','cal_*.npz'),recursive=True);dat=[]
    for fp in files:
        z=np.load(fp);dat.append((float(z['rho']),int(z['seed']),z['y'],z['alra_z'],z['neigh']))
    if len(dat)!=30:raise RuntimeError(f'expected 30 calibration files, got {len(dat)}')
    out={'pass_rule':PASS_RULE,'models':{}};folds=[(set(range(10001,10006)),set(range(10006,10011))),(set(range(10006,10011)),set(range(10001,10006)))]
    for train,test in folds:
        Xh=[];Xf=[];Y=[]
        for rho,seed,y,a,n in dat:
            if seed in train:Xh.append(a[:,None]);Xf.append(np.c_[a,n]);Y.append(y)
        yy=np.concatenate(Y);mh=LogisticRegression(C=1000,max_iter=1000).fit(np.vstack(Xh),yy);mf=LogisticRegression(C=1000,max_iter=1000).fit(np.vstack(Xf),yy);key='seeds_'+str(min(test))+'_'+str(max(test));out['models'][key]={'test_seeds':sorted(test),'h_intercept':float(mh.intercept_[0]),'h_coef':[float(x) for x in mh.coef_[0]],'f_intercept':float(mf.intercept_[0]),'f_coef':[float(x) for x in mf.coef_[0]]}
    json.dump(out,open(args.output,'w'),indent=2)
def bin_det(x):
    bins=np.array([0,.001,.01,.05,.10,.20,.50,1.000001]);labs=np.array(['<0.1%','0.1-1%','1-5%','5-10%','10-20%','20-50%','>=50%']);return labs[np.clip(np.digitize(x,bins)-1,0,len(labs)-1)]
def evaluate(args):
    f=event_features(args.input,args.rho,args.seed,args.neg_cap);co=json.load(open(args.coef));model=next(v for v in co['models'].values() if args.seed in v['test_seeds']);y=f['y'];hbf=model['h_intercept']+model['h_coef'][0]*f['alra_z'];fbf=model['f_intercept']+model['f_coef'][0]*f['alra_z']+model['f_coef'][1]*f['neigh'];hs=f['hier_d']+hbf;fs=f['hier_a']+fbf+f['logq'];methods={'ALRA':(f['alra_z'],f['alra_call']),'scGACL':(lg(f['sc_score']),f['sc_call']),'H_ALRA_v2':(hs,hs>0),'full_H_ALRA_v2':(fs,fs>0)};rows=[];thresholds={}
    for m,(s,c) in methods.items():q=calc_metrics(y,s,c);thresholds[m]=(q['th_fpr1'],q['th_fpr5']);q.update(method=m,rho=args.rho,seed=args.seed,n_pos=int((y==1).sum()),n_neg=int((y==0).sum()));rows.append(q)
    pd.DataFrame(rows).to_csv(args.output,index=False);sr=[];db=bin_det(f['det'][f['g']])
    for b in np.unique(db):
        ix=db==b;yy=y[ix]
        for m,(s,c) in methods.items():
            row={'method':m,'rho':args.rho,'seed':args.seed,'stratum_type':'gene_detection','stratum':b,'n':int(ix.sum()),'n_pos':int((yy==1).sum()),'n_neg':int((yy==0).sum())};row['auroc']=roc_auc_score(yy,s[ix]) if len(np.unique(yy))==2 else np.nan;row['auprc']=average_precision_score(yy,s[ix]) if len(np.unique(yy))==2 else np.nan
            if np.any(yy==1):row['tpr_at_global_fpr1']=float(np.mean(s[ix][yy==1]>thresholds[m][0]));row['tpr_at_global_fpr5']=float(np.mean(s[ix][yy==1]>thresholds[m][1]))
            sr.append(row)
    pd.DataFrame(sr).to_csv(args.output+'.strata.csv',index=False);json.dump(f['diag'],open(args.output+'.diag.json','w'),indent=2)
def aggregate(args):
    files=glob.glob(os.path.join(args.metrics_dir,'**','metrics_*.csv'),recursive=True);df=pd.concat([pd.read_csv(x) for x in files if '.strata.' not in x],ignore_index=True)
    if len(df)!=120:raise RuntimeError(f'expected 120 method rows, got {len(df)}')
    df.to_csv(os.path.join(args.output_dir,'metrics_by_mask.csv'),index=False);df.groupby('method')[['auroc','auprc','tpr_fpr1','tpr_fpr5','tpr_default','fpr_default','precision_default']].agg(['mean','std']).to_csv(os.path.join(args.output_dir,'summary.csv'));wide=df.pivot(index=['rho','seed'],columns='method');gains=[]
    for metric in ['auroc','auprc','tpr_fpr1','tpr_fpr5']:
        base=np.maximum(wide[(metric,'ALRA')],wide[(metric,'scGACL')])
        for m in ['H_ALRA_v2','full_H_ALRA_v2']:
            d=wide[(metric,m)]-base;gains.append({'method':m,'metric':metric,'mean_delta':float(d.mean()),'median_delta':float(d.median()),'wins':int((d>0).sum()),'n':len(d),'mean_best_baseline':float(base.mean()),'mean_method':float(wide[(metric,m)].mean()),'residual_error_reduction':float(1-(1-wide[(metric,m)].mean())/(1-base.mean())) if metric in ['auroc','auprc'] and base.mean()<1 else np.nan})
    gd=pd.DataFrame(gains);gd.to_csv(os.path.join(args.output_dir,'paired_gains.csv'),index=False);verdict=[]
    for m in ['H_ALRA_v2','full_H_ALRA_v2']:
        q=gd[gd.method==m].set_index('metric');mm=df[df.method==m];aa=df[df.method=='ALRA'];rho_ok=True
        for rho in RHO_LIST:
            w=df[df.rho==rho].pivot(index='seed',columns='method',values=['tpr_fpr1','auprc']);rho_ok &= bool(((w[('tpr_fpr1',m)]-np.maximum(w[('tpr_fpr1','ALRA')],w[('tpr_fpr1','scGACL')])).mean()>0) and ((w[('auprc',m)]-np.maximum(w[('auprc','ALRA')],w[('auprc','scGACL')])).mean()>0))
        crit={'tpr1':bool(q.loc['tpr_fpr1','mean_delta']>=PASS_RULE['tpr_fpr1_min_delta'] and q.loc['tpr_fpr1','wins']>=PASS_RULE['paired_wins_min']),'tpr5':bool(q.loc['tpr_fpr5','mean_delta']>=PASS_RULE['tpr_fpr5_min_delta'] and q.loc['tpr_fpr5','wins']>=PASS_RULE['paired_wins_min']),'auprc':bool(q.loc['auprc','residual_error_reduction']>=PASS_RULE['auprc_residual_reduction_min'] and q.loc['auprc','wins']>=PASS_RULE['paired_wins_min']),'auroc':bool(q.loc['auroc','residual_error_reduction']>=PASS_RULE['auroc_residual_reduction_min'] and q.loc['auroc','wins']>=PASS_RULE['paired_wins_min']),'default':bool(mm.tpr_default.mean()>=aa.tpr_default.mean() and mm.fpr_default.mean()<=aa.fpr_default.mean()+PASS_RULE['default_fpr_margin_vs_alra']),'rho_consistency':bool(rho_ok)};verdict.append({'method':m,'PASS':all(crit.values()),**crit})
    vd=pd.DataFrame(verdict);vd.to_csv(os.path.join(args.output_dir,'verdict.csv'),index=False);strata=glob.glob(os.path.join(args.metrics_dir,'**','*.strata.csv'),recursive=True)
    if strata:pd.concat([pd.read_csv(x) for x in strata],ignore_index=True).to_csv(os.path.join(args.output_dir,'stratified_metrics.csv'),index=False)
    diags=glob.glob(os.path.join(args.metrics_dir,'**','*.diag.json'),recursive=True);diag=[json.load(open(x)) for x in diags];json.dump(diag,open(os.path.join(args.output_dir,'diagnostics.json'),'w'),indent=2);prov=diag[0]['provenance'] if diag else {};lines=['# Full-gene H-ALRA detector benchmark','',f"- Gene universe: **{prov.get('n_genes_all','NA')} genes** (no pre-benchmark gene filtering).",f"- Cells: **{prov.get('n_cells','NA')}**.",'- 30 paired masks: 3 retention levels × 10 seeds.','- All methods receive the same full gene universe and same masks.','- scGACL >95%-zero invalidity is counted as method coverage failure; genes are not removed from the common universe.','- Strict lineage-restricted negatives use a deterministic cap of 100 zero coordinates per gene×lineage; every qualifying gene remains represented.','','## Preregistered PASS/FAIL rule','```json',json.dumps(PASS_RULE,indent=2),'```','','## Mean metrics','',df.groupby('method')[['auroc','auprc','tpr_fpr1','tpr_fpr5','tpr_default','fpr_default']].mean().to_markdown(),'','## Paired gains versus better ALRA/scGACL baseline','',gd.to_markdown(index=False),'','## Verdict','',vd.to_markdown(index=False),'','## Provenance','```json',json.dumps(prov,indent=2),'```'];open(os.path.join(args.output_dir,'REPORT.md'),'w').write('\n'.join(lines))
def main():
    p=argparse.ArgumentParser();sub=p.add_subparsers(dest='cmd',required=True);a=sub.add_parser('features');a.add_argument('--input',required=True);a.add_argument('--rho',type=float,required=True);a.add_argument('--seed',type=int,required=True);a.add_argument('--output',required=True);a.add_argument('--neg-cap',type=int,default=100);a.add_argument('--cal-per-class',type=int,default=10000);a=sub.add_parser('calibrate');a.add_argument('--cal-dir',required=True);a.add_argument('--output',required=True);a=sub.add_parser('evaluate');a.add_argument('--input',required=True);a.add_argument('--rho',type=float,required=True);a.add_argument('--seed',type=int,required=True);a.add_argument('--coef',required=True);a.add_argument('--output',required=True);a.add_argument('--neg-cap',type=int,default=100);a=sub.add_parser('aggregate');a.add_argument('--metrics-dir',required=True);a.add_argument('--output-dir',required=True);z=p.parse_args();os.makedirs(os.path.dirname(z.output) if hasattr(z,'output') and os.path.dirname(z.output) else '.',exist_ok=True);{'features':save_feature_mode,'calibrate':calibrate,'evaluate':evaluate,'aggregate':aggregate}[z.cmd](z)
if __name__=='__main__':main()
