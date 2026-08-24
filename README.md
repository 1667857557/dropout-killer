# DropoutKiller

`DropoutKiller` is an R package for **selective** scRNA-seq dropout recovery. It does not impute an entire expression matrix. It first identifies observed zeros with strong membership-local evidence of technical dropout, then modifies only those events.

Recovery is expression-derived only: **no PPI, pathway, GRN, or other external biological prior is used**.

## Workflow

```text
normalized scRNA-seq expression X (genes x cells)
                 |
                 +-----------------------------+
                 |                             |
                 v                             v
        low-dimensional embedding       broad biological strata
                 |                       (e.g. major cell type)
                 +-------------+---------------+
                               |
                               v
              SuperCell-style membership construction
       kNN graph -> walktrap hierarchy -> gamma-controlled cut
                               |
                               v
                     membership-local ALRA
                               |
                               v
          zero gate + negative-tail confidence score
                               |
                    only observed X[g,c] = 0
                               |
                               v
                   high-confidence dropout mask
                               |
                               v
        membership-local Gaussian neighbor borrowing
             distance-derived weights, normalized
                               |
                               v
                    selective replacement only
                               |
                               v
                     recovered expression
```

## Recovery model

For a masked event `(g,c)`, donors are restricted to the same membership. For donor cell `j`,

```text
d_cj^2 = ||z_c - z_j||_2^2
w_cj ∝ exp(-d_cj^2 / sigma_c^2)
```

By default `sigma_c` is estimated as the median positive distance among the selected neighbors. With `neighbor_positive_only = TRUE`, only cells with observed `X[g,j] > 0` contribute and the weights are renormalized over those donors:

```text
Xhat[g,c] = sum_j w_cj X[g,j] I(X[g,j] > 0) /
            sum_j w_cj I(X[g,j] > 0)
```

There is no fixed cell/prior mixing coefficient and no biological-prior branch.

## Installation

```r
remotes::install_github("1667857557/dropout-killer", upgrade = "never")
```

## Matrix workflow

```r
library(DropoutKiller)

fit <- dropout_killer(
  x = x,
  embedding = pca,
  group = major_cell_type,
  split_by = condition,
  gamma = 20,
  k_knn = 5,
  rank = "auto",
  quantile_prob = 0.001,
  threshold = 0.95,
  neighbor_k = 30,
  neighbor_positive_only = TRUE
)

fit
fit$expression
fit$events
validate_dropout_result(fit, x)
```

If condition-specific states must never borrow from one another, pass condition to `split_by`. If cross-condition borrowing within the same trusted biological state is intended, leave `split_by = NULL`.

## Direct recovery for a supplied mask

```r
recovered <- recover_dropout_expression(
  x = x,
  mask = dropout_mask,
  membership = membership,
  embedding = pca,
  neighbor_k = 30,
  neighbor_positive_only = TRUE
)
```

## Seurat workflow

```r
obj <- dropout_killer_seurat(
  object = obj,
  assay = "RNA",
  slot = "data",
  reduction = "pca",
  dims = 1:20,
  group_by = "major_cell_type",
  split_by = "condition",
  new_assay = "DropoutKiller"
)
```

Recovered values are continuous and are stored as **assay data**, not as raw counts.

## Main API

- `build_supercell_membership()`: graph-based membership construction.
- `local_alra_detect()`: membership-local low-rank zero detection.
- `select_dropout_mask()`: high-confidence zero-only sparse mask.
- `weighted_neighbor_prediction()`: membership-local Gaussian neighbor prediction.
- `recover_dropout_expression()`: selective neighbor-only recovery for a supplied mask.
- `dropout_killer()`: end-to-end workflow.
- `dropout_killer_seurat()`: Seurat wrapper.
- `validate_dropout_result()`: validate zero-preserving invariants.

See `inst/ALGORITHM.md` for the mathematical contract.

## Important trade-offs

- Hard `group`/`split_by` boundaries should use trusted broad labels. Over-fragmentation reduces membership size and local low-rank power.
- `threshold = 0.95` is a null-tail evidence threshold, not an FDR or calibrated posterior probability.
- Positive-only neighbor borrowing estimates a conditional-positive local mean and can be upward-biased. Set `neighbor_positive_only = FALSE` for the literal all-neighbor Gaussian average.
- Recovery is one-pass: recovered values are never fed back into membership construction or dropout detection.

## Methodological basis

The membership engine adapts the SuperCell graph/coarse-graining design: kNN graph, walktrap hierarchy, and gamma-controlled coarse graining. The zero gate and automatic-rank logic adapt ALRA's low-rank and negative-tail ideas. The recovery stage is deliberately limited to membership-local, distance-weighted cell borrowing.
