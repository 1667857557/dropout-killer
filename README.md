# DropoutKiller

`DropoutKiller` is an R package for **selective** scRNA-seq dropout recovery. It does not impute an entire matrix. It first asks whether an observed zero has strong local evidence of being a technical dropout and modifies only events that pass that gate.

The implementation follows the framework in this repository while tightening several statistical and engineering points:

1. **SuperCell-style membership**: within user-supplied biological strata (typically major cell type), cells are connected by a Euclidean kNN graph in a low-dimensional embedding; walktrap hierarchical clustering is cut to approximately `n/gamma` memberships. Large strata can use the same anchor/centroid approximation idea as SuperCell.
2. **Local ALRA gate**: each membership is decomposed independently. For gene `g`, the low-rank reconstruction is compared with an adaptive threshold derived from its negative tail. Only *observed zeros* above this threshold are candidates.
3. **Dropout confidence**: negative low-rank values estimate a symmetric biological-zero error distribution. The reported confidence is one-sided evidence against that null; it is deliberately **not described as a calibrated posterior dropout probability**.
4. **Cell-space recovery**: candidate expression is borrowed only from cells in the same membership, with Gaussian weights in the embedding and, by default, only from neighbors where that gene is observed (`>0`).
5. **Gene-space prior**: optional PPI/pathway/GRN adjacency matrices are row-normalized and applied in membership-local standardized expression space to avoid directly averaging genes with incompatible expression scales.
6. **Selective replacement**: when both components are available the default is `75%` cell borrowing + `25%` gene prior. If one component is unavailable, weights are renormalized over the available component. Observed non-zero values are never overwritten.

## Installation

```r
remotes::install_github("1667857557/dropout-killer")
```

## Core matrix workflow

`x` should be a non-negative normalized expression matrix with genes in rows and cells in columns. `embedding` should have cells in rows (typically PCA).

```r
library(DropoutKiller)

fit <- dropout_killer(
  x = x,
  embedding = pca,
  group = major_cell_type,
  split_by = condition,       # optional hard boundary
  gamma = 20,
  k_knn = 5,
  rank = "auto",
  quantile_prob = 0.001,
  threshold = 0.95,
  alpha = 0.75,
  neighbor_k = 30,
  neighbor_positive_only = TRUE,
  gene_networks = list(ppi, pathway, grn)
)

fit
fit$expression
fit$events
validate_dropout_result(fit, x)
```

If condition-specific differences must not be borrowed across conditions, pass the condition vector to `split_by`. If the purpose is to allow local cross-condition borrowing within the same annotated cell state, leave `split_by = NULL` and interpret downstream condition comparisons accordingly.

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

The recovered values are continuous and are stored as **assay data**, not raw counts. For count-based models, retain the original counts and use the recovery output only where a continuous expression representation is appropriate.

## Why the membership implementation is not a literal copy of SuperCell

SuperCell builds a kNN graph, clusters it, and can split metacells by annotation after graph construction. DropoutKiller uses the same graph/coarse-graining logic but applies supplied biological strata **before** graph construction. This prevents cross-stratum graph edges from influencing membership boundaries, which is a stricter requirement for selective expression borrowing.

## Main objects

- `build_supercell_membership()`: graph-based membership construction.
- `local_alra_detect()`: membership-local low-rank detection with adaptive negative-tail thresholding.
- `select_dropout_mask()`: high-confidence sparse mask.
- `weighted_neighbor_prediction()`: masked, membership-constrained cell prediction.
- `gene_network_from_edges()` / `prepare_gene_network()` / `combine_gene_prior()`: edge-list construction and aligned PPI/pathway/GRN priors.
- `gene_prior_prediction()`: masked gene-space prediction.
- `recover_dropout_expression()`: selective hybrid recovery.
- `dropout_killer()`: end-to-end pipeline.
- `dropout_killer_seurat()`: Seurat wrapper.

See `inst/ALGORITHM.md` for the mathematical contract and implementation invariants.

## Important design trade-offs

- Hard `group`/`split_by` boundaries should use trusted broad labels. Noisy or over-fragmented annotations reduce membership size and therefore local low-rank power.
- The 0.95 confidence gate is a null-tail evidence threshold, not a calibrated error rate. It should be validated by synthetic masking on the target dataset when strong quantitative claims depend on recovered values.
- Positive-only neighbor borrowing avoids dilution by unresolved zeros but estimates expression conditional on observed donors and can be upward-biased. Set `neighbor_positive_only = FALSE` to use all local neighbors.
- PPI/pathway/GRN information is an optional prior, not evidence that a target is expressed. It is used only after the zero event passes the expression-derived dropout gate.
- Recovery is one-pass: recovered values are never fed back into membership construction or dropout detection, avoiding iterative self-reinforcement.


## Methodological references

The membership engine is a critical adaptation of the SuperCell graph/coarse-graining design (kNN graph, walktrap hierarchy, gamma-controlled cut). The zero gate and automatic-rank logic are based on ALRA's low-rank/negative-tail construction. DropoutKiller does not reproduce either package verbatim; the implementation contracts above define where behavior intentionally differs for selective recovery.
