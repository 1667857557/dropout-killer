#!/usr/bin/env python3
"""Calibration-only fast path for the strict all-gene H-ALRA benchmark.

The original `features` command calls the complete event_features() pipeline,
including scGACL EM and the hierarchical/NB observation model, even though OOF
calibration uses only ALRA z and frozen-neighbor evidence.  This script computes
exactly those two calibration features on the exact same masks/negative controls.
It changes scheduling only, not any fitted coefficient, detector score, or final
evaluation formula.
"""
import argparse
import json
import os
import time
import numpy as np

# Importing the strict wrapper installs the strict all-positive mask generator and
# exact scGACL invalidity behavior into the audited base module.  This fast path
# uses only the strict mask patch; final evaluation still runs the strict wrapper.
import full_gene_h_alra_benchmark_strict as strict

base = strict.base


def calibration_features(inp, rho, seed, neg_cap=100, cal_per_class=10000):
    t0 = time.time()
    C0, genes, meta, lab, cli, emb, mem = base.read_input(inp)
    prov, det = base.provenance(C0, lab)
    if not prov['all_pass']:
        raise RuntimeError('PROVENANCE FAILURE: ' + json.dumps(prov))

    npos0 = base.supports(C0, cli)
    ng, nc, nl = base.build_neg_controls(C0, cli, npos0, cap=neg_cap)
    C, mg, mc, mv = base.make_mask(C0, rho, seed, cli, npos0)

    lost = np.bincount(mc, weights=mv, minlength=C0.shape[1])
    lib = np.maximum(np.asarray(C0.sum(0)).ravel() - lost, 1)
    eg = np.r_[mg, ng].astype(np.int32)
    ec = np.r_[mc, nc].astype(np.int32)
    y = np.r_[np.ones(len(mg), np.int8), np.zeros(len(ng), np.int8)]

    alra_z, _alra_call, rank = base.alra_events(C, lib, eg, ec, seed)
    W = base.build_W(emb, cli, 30)
    neigh = base.neighbor_events(C, W, eg, ec)

    # Byte-for-byte equivalent sampling rule to base.sample_cal().
    rng = np.random.default_rng(seed + 20260907)
    pp = np.where(y == 1)[0]
    nn = np.where(y == 0)[0]
    k = min(int(cal_per_class), len(pp), len(nn))
    ix = np.r_[rng.choice(pp, k, False), rng.choice(nn, k, False)]

    diag = {
        'rho': float(rho),
        'seed': int(seed),
        'calibration_only_fast_path': True,
        'n_genes_all': int(C0.shape[0]),
        'n_cells': int(C0.shape[1]),
        'masked_n': int(len(mg)),
        'negative_n': int(len(ng)),
        'calibration_n_per_class': int(k),
        'alra_rank_auto': int(rank),
        'baseline_zero_fraction_all': float(1 - C0.nnz / (C0.shape[0] * C0.shape[1])),
        'masked_zero_fraction_all': float(1 - C.nnz / (C.shape[0] * C.shape[1])),
        'runtime_sec': float(time.time() - t0),
        'provenance': prov,
        'neg_cap_per_gene_lineage': int(neg_cap),
        'skipped_as_mathematically_irrelevant_to_calibration': [
            'scGACL EM', 'hierarchical prior', 'NB/Poisson zero likelihood'
        ],
    }
    return y[ix], alra_z[ix], neigh[ix], diag


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', required=True)
    p.add_argument('--rho', type=float, required=True)
    p.add_argument('--seed', type=int, required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--neg-cap', type=int, default=100)
    p.add_argument('--cal-per-class', type=int, default=10000)
    args = p.parse_args()
    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    y, a, n, diag = calibration_features(
        args.input, args.rho, args.seed, args.neg_cap, args.cal_per_class
    )
    np.savez_compressed(
        args.output,
        y=y,
        alra_z=a,
        neigh=n,
        rho=float(args.rho),
        seed=int(args.seed),
    )
    with open(args.output + '.diag.json', 'w') as fh:
        json.dump(diag, fh, indent=2)


if __name__ == '__main__':
    main()
