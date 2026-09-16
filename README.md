# DropoutKiller

DropoutKiller detects **recoverable RNA zero events** and selectively reconstructs only those calls. The current production framework is:

```text
raw RNA counts
  -> library-size normalization + log1p
  -> SuperCell hierarchy geometry within broad biological classes
  -> Lean recoverable-zero detector
       [ALRA low-rank evidence + local hierarchy support + membership context]
  -> calibrated sparse zero mask
  -> target-safe P1 stabilized-state recovery
  -> recovered normalized expression
```

Observed nonzero expression is never overwritten. The package does **not** use scGACL, a Gamma-Normal mixture, or a gene-by-cell-type dropout posterior.

## Primary usage

```r
# broad_type must represent broad biological lineages, not fine subclusters.
out <- dropout_killer_seurat(
  object,
  group_by = "broad_type",
  return_result = TRUE
)

out$result$settings$detection_method
model <- out$result$detection$model

# Reuse an independently validated detector model of the same modality.
out2 <- dropout_killer_seurat(
  object,
  group_by = "broad_type",
  lean_model = model,
  return_result = TRUE
)
```

For a gene-by-cell RNA count matrix, use:

```r
res <- dropout_killer(
  counts,
  embedding = recovery_embedding,
  group = broad_group
)
```

The supplied matrix embedding controls recovery geometry. RNA-only detection rebuilds its own normalized-RNA SVD geometry so calibration and inference use the same feature-generating pipeline.

## 1. Working expression scale

Raw RNA counts are normalized to 10,000 counts per cell and transformed with `log1p`:

\[
Y_{gc}=\log\left(1+10^4\frac{X_{gc}}{\sum_h X_{hc}}\right).
\]

Recovery produces continuous values on this normalized log-expression scale. The Seurat wrapper therefore writes the recovered matrix to the `data` slot of a new assay rather than presenting it as integer counts.

## 2. Broad biological boundaries and hierarchy

`group` / `group_by` defines hard broad biological boundaries. Examples are major lineages such as tumor, T/NK, B/plasma, myeloid, fibroblast, or endothelial cells. Fine cell-state labels should not be supplied as hard groups.

Within each broad class, DropoutKiller builds a kNN graph, retains the full Walktrap merge history, and cuts the graph at the SuperCell-style target

\[
K=\max\left(1,\left\lfloor\frac{n}{\gamma}\right\rfloor\right),
\]

with `gamma = 150` by default, subject to connected-component constraints. No cross-class borrowing is introduced.

For a query cell, support cells are taken from its own final membership first. Only when that support is insufficient are cells borrowed from the closest connected memberships in the retained hierarchy.

## 3. Lean recoverable-zero detector

For every observed RNA zero, the production detector combines three cell-specific features.

### ALRA low-rank evidence

The normalized RNA matrix is reconstructed with a data-driven low-rank SVD. The detector uses a standardized gene-specific margin `z` derived from the reconstructed value and the negative reconstruction scale. This is latent-expression evidence, not a dropout posterior.

### Local hierarchy support

Up to `neighbor_k = 30` support cells are selected from the query cell's own membership and, only on shortage, nearby memberships in the retained hierarchy. The weighted fraction of positive neighbors is transformed with the benchmark smoothing rule

\[
L=\operatorname{logit}\left(\frac{s+0.01}{1.02}\right).
\]

### Membership context

For gene \(g\), query membership \(m\), and its broad class \(G\), the context term is

\[
H=\frac{n_m}{n_m+50}
\left[
\operatorname{logit}\left(\frac{k_m+0.5}{n_m+1}\right)-
\operatorname{logit}\left(\frac{k_G+0.5}{n_G+1}\right)
\right].
\]

This asks whether the gene is locally active in the query state relative to its broader lineage while shrinking small memberships toward zero context evidence.

The three features are combined by a class-balanced logistic model:

\[
\eta=\beta_0+\beta_1z+\beta_2L+\beta_3H.
\]

A zero is called only when its **linear score** is strictly above the model's calibrated threshold. `plogis(eta)` is stored for convenience but is not interpreted as a calibrated biological dropout probability.

## 4. Calibration

When `lean_model` is absent, raw UMI counts are calibrated by artificial thinning. Defaults are `q = 0.5` UMI retention and seeds `10001:10003`. Each thinned matrix rebuilds the same detector geometry used for inference.

Positive training events are original positive counts that become zero after thinning. Proxy negatives are original zeros whose gene prevalence is at most 0.005 in the current broad class and at least 0.20 in another broad class. Calibration requires at least two broad classes and stops rather than silently relaxing these constraints.

The default `fpr = 0.01` is a **training proxy-FPR target**. It is not a biological false-discovery rate and must not be reported as one. Independent masks or datasets are required for performance evaluation.

## 5. RNA + ATAC multiome geometry

For an unambiguous paired RNA+ATAC Seurat object, `dropout_killer_seurat()` automatically uses `Supercell_hierarchy_Lean_membership_WNN` unless `modality = "rna"` explicitly opts out.

Every WNN build recomputes:

- RNA normalization, variable features, scaling, and PCA;
- ATAC TF-IDF, top features, and LSI;
- PCA dimensions 1:40 and LSI dimensions 2:40, reduced as needed for small inputs;
- `FindMultiModalNeighbors()` separately inside each broad class.

LSI component 1 is always excluded. The hierarchy graph follows the SuperCell 2.0 multimodal-neighbor kernel

\[
A_{ij}=1-2d_{ij}^2,
\]

followed by symmetrization. It does not substitute the Seurat `wsnn` shared-neighbor graph or concatenated PCA/LSI coordinates.

ATAC therefore contributes to **biological neighborhood geometry**. It is not directly used to synthesize RNA expression values.

## 6. Selective P1 stabilized-state recovery

Detection and recovery are separate estimands. Only coordinates in the sparse detector mask are eligible for reconstruction.

The default recovery engine is `p1_stabilized_state`:

1. target genes are assigned to deterministic cross-fitting folds;
2. genes in the current target fold are excluded from predictor-state construction;
3. non-target standardized expression is smoothed once over the retained hierarchy/embedding geometry;
4. a low-dimensional cell-state representation is learned from those target-safe predictors;
5. positive target-gene donors fit a ridge state model;
6. analytic leave-one-out diagnostics determine shrinkage and bias calibration;
7. unsupported target states fall back conservatively rather than forcing a factor prediction.

The default production settings are five target-gene folds, ridge penalty 2, support-adaptive factor rank, bias shrinkage `bias_kappa = 10`, and predictor smoothing 0.25.

`recover_dropout_expression()` exposes the P1 engine directly for a supplied zero mask. Older recovery estimators remain internal comparison components used for ablation; their former standalone wrapper APIs are no longer part of the package interface.

## 7. Public API

The supported public surface is intentionally small:

- `dropout_killer()` — matrix workflow;
- `dropout_killer_seurat()` — Seurat RNA or RNA+ATAC workflow;
- `fit_lean_detector()` and `calibrate_lean_detector()` — detector model fitting/calibration;
- `supercell_lean_detect()` — detection-only execution;
- `build_supercell_membership()` and `build_wnn_supercell()` — inspectable geometry builders;
- `recover_dropout_expression()` — selective recovery for a supplied mask;
- `sample_dropout_expression()` — uncertainty-aware completed-matrix draws when predictive variance is available;
- `membership_summary()` and `validate_dropout_result()` — diagnostics.

Historical scGACL, local/global ALRA detector routes, compatibility aliases, and standalone comparator wrappers are not part of the current API.

## 8. Validation boundary

A stronger post-recovery correlation is not evidence that recovery is correct when correlation information contributed to prediction. Evaluation should use independent count thinning or held-out positive events and report detection and recovery separately. In particular, compare recoverable-zero precision/recall, conditional reconstruction error, downstream differential-expression stability, clustering/marker preservation, and uncertainty coverage on masks that were not used to fit the detector threshold.
