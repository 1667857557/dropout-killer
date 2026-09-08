#!/usr/bin/env python3
"""Strict all-gene wrapper for the full H-ALRA detector benchmark.

This wrapper deliberately removes every pre-benchmark gene/support gate from the
synthetic dropout generator. Every originally non-zero UMI coordinate in the
36,601-gene RNA matrix is eligible for count-dependent thinning:

    P(k -> 0) = (1 - rho) ** k

Genes with no observed non-zero UMI remain in the common 36,601-gene method
universe but cannot, by definition, yield a synthetic positive coordinate.
Native zeros are never relabelled wholesale as biological negatives.
"""
import numpy as np
from scipy import sparse
import full_gene_h_alra_benchmark as base


def strict_make_mask(C0, rho, seed, cli=None, npos0=None):
    """Thin all observed positive coordinates without gene/support/cell caps."""
    rng = np.random.default_rng(seed)
    coo = C0.tocoo(copy=False)
    g = coo.row.astype(np.int32, copy=False)
    c = coo.col.astype(np.int32, copy=False)
    v = coo.data.astype(np.int32, copy=False)

    # No eligibility gate: every observed positive UMI coordinate is eligible.
    p = np.power(1.0 - float(rho), v.astype(np.float64))
    sel = rng.random(len(v)) < p

    mg = g[sel]
    mc = c[sel]
    mv = v[sel]
    keep = ~sel
    C = sparse.csc_matrix(
        (v[keep], (g[keep], c[keep])), shape=C0.shape, dtype=np.int32
    )
    return C, mg, mc, mv


def exact_scgacl_params(C, cli):
    """Mathematically equivalent scGACL fit with exact invalid-group pruning.

    The authors' detector declares a gene/group invalid before EM when its
    estimated zero rate is >0.95 (or it has no positives). We therefore avoid
    calling EM for those groups, but retain them in the common gene universe;
    downstream base code scores invalid scGACL entries as its coverage failure.
    For every group that can be valid under this gate, the unchanged authors'
    translated EM implementation (base.sc_fit_sparse) is called.
    """
    G, N = C.shape
    L = 6
    d = np.full((G, L), np.nan, np.float32)
    lam = np.full((G, L), np.nan, np.float32)
    invalid = np.ones((G, L), dtype=bool)

    for li in range(L):
        ids = np.where(cli == li)[0]
        R = C[:, ids].tocsr()
        n = len(ids)
        npos = np.diff(R.indptr)
        # rate = (n-npos)/n; original invalid gate is rate > 0.95.
        candidates = np.flatnonzero((npos > 0) & (((n - npos) / n) <= 0.95))
        for g in candidates:
            vals = R.data[R.indptr[g]:R.indptr[g + 1]]
            dd, ll, ii = base.sc_fit_sparse(vals, n)
            d[g, li] = dd
            lam[g, li] = ll
            invalid[g, li] = ii

    R = C.tocsr()
    dg = np.full(G, np.nan, np.float32)
    glam = np.full(G, np.nan, np.float32)
    invg = np.ones(G, dtype=bool)
    npos = np.diff(R.indptr)
    candidates = np.flatnonzero((npos > 0) & (((N - npos) / N) <= 0.95))
    for g in candidates:
        vals = R.data[R.indptr[g]:R.indptr[g + 1]]
        dd, ll, ii = base.sc_fit_sparse(vals, N)
        dg[g] = dd
        glam[g] = ll
        invg[g] = ii
    return d, lam, dg, glam, invalid, invg


_original_event_features = base.event_features


def strict_event_features(inp, rho, seed, neg_cap=100):
    f = _original_event_features(inp, rho, seed, neg_cap)
    det = np.asarray(f['det'])
    pos = f['y'] == 1
    masked_genes = np.unique(f['g'][pos])
    n_observed = int(np.sum(det > 0))
    G = int(len(det))
    f['diag'].update({
        'benchmark_gene_universe': 'all_RNA_features',
        'strict_all_gene_masking': True,
        'mask_formula': 'P(k_to_0)=(1-rho)^k for every originally nonzero coordinate',
        'gene_support_filter': None,
        'per_cell_mask_cap': None,
        'per_gene_lineage_retention_floor': None,
        'global_mask_cap': None,
        'n_genes_total': G,
        'n_genes_excluded_by_benchmark_filter': 0,
        'n_genes_with_any_observed_positive': n_observed,
        'n_genes_all_observed_zero': G - n_observed,
        'n_genes_with_masked_positive_this_mask': int(masked_genes.size),
        'fraction_observed_positive_genes_masked_this_mask': (
            float(masked_genes.size / n_observed) if n_observed else np.nan
        ),
    })
    if G != 36601:
        raise RuntimeError(f'STRICT ALL-GENE FAILURE: expected 36601 genes, got {G}')
    if f['diag']['n_genes_excluded_by_benchmark_filter'] != 0:
        raise RuntimeError('STRICT ALL-GENE FAILURE: a benchmark gene filter was applied')
    return f


# Monkey-patch only the two operations that need correction/optimization.
# All ALRA, hierarchy, neighbor, OOF calibration, metrics and PASS/FAIL code
# remain exactly the audited base benchmark implementation.
base.make_mask = strict_make_mask
base.scgacl_params = exact_scgacl_params
base.event_features = strict_event_features

if __name__ == '__main__':
    base.main()
