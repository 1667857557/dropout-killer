# DropoutKiller

`DropoutKiller` is an R package for **selective** scRNA-seq dropout recovery. It never overwrites observed non-dropout coordinates and does not use external PPI/pathway/GRN priors.

Version 0.6 keeps detection and recovery as separate statistical problems. Following the real-data artificial-dropout benchmark, production zero detection is now **native-ALRA-style global detection within each major cell class**, while SuperCell memberships and hierarchy remain local recovery structures.

```text
raw expression X + major cell class + biological embedding
                         |
                         v
          ALRA library-size + log normalization
                         |
                         v
      native ALRA rank-k reconstruction per cell class
                         |
                         v
        gene-wise |Q_0.001(low-rank)| zero gate
                         |
                         v
                 selective dropout mask
                         |
          +--------------+----------------+
          |                               |
          v                               v
  SuperCell gamma cut              retained tree
          |                               |
          +---------------+---------------+
                          |
                          v
                 local recovery engines
```

## Detection

The default detector is now:

```r
detection_method = "alra_global_by_group"
```

`group` defines the major cell-class detection blocks. Each class is processed independently, but final SuperCell memberships do **not** fragment zero detection. If `group = NULL`, all cells form one global ALRA block.

For a normalized cell-by-gene block `A`, automatic rank selection follows the original ALRA singular-value-spacing heuristic. With up to `K = 100` singular values, the spacing sequence is

```text
d_i = sigma_i - sigma_{i+1}
```

and rank is selected when a spacing exceeds the noise-tail mean by `rank_z = 6` standard deviations. For cell classes with fewer than 100 available singular values, `K` and the noise-tail start are reduced to the available dimensions while retaining at least five tail spacings when possible.

The final randomized rank-k reconstruction is thresholded gene-wise using the original ALRA rule:

```text
tau_g = |Q_0.001(Ahat_.g)|
call_gc = (A_gc == 0) & (Ahat_gc > tau_g)
```

There is **no second `confidence >= 0.95` gate** for this detector. `threshold` remains relevant only to the historical confidence-based local detector.

The previous membership-local engines remain available explicitly:

```r
detection_method = "eb_zero_null"
detection_method = "alra_quantile"
```

The benchmark motivating the default change showed that strict membership-local detection could lose substantial sensitivity in 8–49-cell memberships, whereas global ALRA retained high artificial-dropout recall. The detector and recovery scopes are therefore deliberately separated.

## Why retain the SuperCell hierarchy?

A SuperCell membership is only a cut through a walktrap hierarchy:

```text
biological embedding -> kNN graph -> walktrap tree -> gamma cut
```

Treating all cells in the final cut as exchangeable discards information already present in that tree. Version 0.6 therefore keeps the `cluster_walktrap()` hierarchy and a compact tree index inside `DropoutKillerMembership`.

For production recovery the final membership is the computational and biological borrowing block:

- **explicit `group` / `split_by`**: additional absolute biological recovery boundary;
- **final SuperCell membership**: absolute recovery borrowing block;
- **walktrap hierarchy + original embedding**: continuous weighting **inside** that membership.

Detection differs intentionally: `group` is the major cell-class ALRA block, while `split_by` and final membership do not fragment the default detector.

## Tree-local donor weighting

For query cell `c` and candidate donor `j` inside the same final membership, the engine uses

```text
w_cj = exp[-alpha d_tree(c,j)/tau
           -(1-alpha) d_embed(c,j)^2/(2 h_c^2)]
```

where:

- `alpha = tree_weight`;
- `d_tree` is a normalized walktrap lowest-common-ancestor distance;
- `d_embed` is distance in the original PCA/Harmony/WNN-like biological embedding;
- `h_c` is an adaptive bandwidth defined by the `local_k`-th candidate neighbor inside the membership.

If hierarchy distance is unavailable for a pair, the pair uses the embedding component instead of inventing a tree distance.

The sparse weight matrix is built once per final membership block and reused across all target genes. `W^2`, row-weight sums, and cell-level tree/embedding distance diagnostics are cached at the same time.

## Query-specific positive baseline

Once a coordinate has already been selected as technical dropout, the positive-expression magnitude baseline becomes query-specific:

```text
mu_gc = sum_j w_cj X_gj / sum_j w_cj
```

using reliable positive donors only.

Its Kish effective donor size is

```text
n_eff,gc = (sum_j w_cj)^2 / sum_j w_cj^2
```

and event diagnostics also retain weighted local positive variance and local positive prevalence.

Unlike a whole-membership positive mean, `mu_gc` varies continuously with biological state inside a membership.

For speed, local means, variances, positive prevalence, and `n_eff` are evaluated in gene batches with matrix multiplication rather than by looping over every dropout event.

## Batched residual coexpression

The tree-local engine does not ask coexpression factors to predict the entire target expression. It first defines leave-one-out local residuals

```text
r_gj = X_gj - mu_gj^(-j)
```

and predicts only what the tree/embedding-local positive baseline cannot explain.

All current target genes are excluded from factor-state learning, preserving the no-target-leakage path

```text
X_cell,-targets -> factor state -> target residual
```

The original implementation fit one weighted ridge for every dropout event, which made runtime scale with millions of high-confidence coordinates. The batched implementation instead fits **one residual model per target gene and final membership**.

Let `Q_g` be the dropout query cells for gene `g`. Donor influence in the single gene-level ridge is still determined by proximity to the actual query cells through the aggregate weight

```text
wbar_gj = sum_{c in Q_g} w_cj
```

and

```text
beta_g = (F' Wbar_g F + P)^-1 F' Wbar_g r_g
```

so donors close to the target dropout cells receive more influence, but the same gene is not refit separately for every event.

Held-out residual predictions estimate a gene-level `q_pred`. Query-specific local information contributes

```text
q_info,gc = n_eff,gc / (n_eff,gc + local_info_kappa)
q_final,gc = q_pred,g * q_info,gc
```

and the deployed prediction is

```text
Xhat_gc = max(0, mu_gc + q_final,gc * rhat_gc)
```

If the factor model has no held-out gain, recovery falls back to the **query-specific local positive mean**, not the whole-membership mean.

## Why the batched implementation is faster

For a membership containing `n_m` cells, biological geometry is constructed once. Recovery then scales approximately with:

```text
one W_m construction
+ batched gene x cell local statistics
+ one small ridge solve per target gene x membership
+ final indexing at selected dropout coordinates
```

rather than:

```text
one small ridge solve per dropout event
```

Thus millions of selected dropout events increase output/indexing work but no longer create millions of independent regression fits.

## Predictive uncertainty and DV

Tree-local fallback variance retains both positive-expression outcome variance and local-mean estimation uncertainty. For a supported residual model, held-out error is evaluated using the same query-specific deployed `q_final` used in the recovered mean.

The working variance combines residual error, local-mean uncertainty, and weighted factor leverage. It remains an approximate predictive variance and should be checked with held-out coverage.

For positive-conditional recovery, `sample_dropout_expression()` uses Gamma moment matching. For predictive mean `m > 0` and variance `v`,

```text
shape = m^2 / v
scale = v / m
X ~ Gamma(shape, scale)
```

so `E[X] = m`, `Var(X) = v`, and draws remain positive.

For differential variability,

```text
Var(lambda_g | Y)
  = Var_i(E[lambda_ig | Y])
  + E_i(Var(lambda_ig | Y))
```

so the deterministic mean matrix alone is not a complete DV representation. Use `prediction_sd`, `predictive_variance`, or repeated completed draws for DV/covariance sensitivity analysis.

## Matrix workflow

```r
library(DropoutKiller)

fit <- dropout_killer(
  x = x,
  embedding = pca,
  group = major_cell_type,
  split_by = condition,
  gamma = 150,
  detection_method = "alra_global_by_group",
  quantile_prob = 0.001,
  recovery_method = "tree_local_factor",
  factor_target = "positive",
  tree_weight = 0.5,
  local_k = 30,
  candidate_k = 100,
  min_effective_donors = 5,
  local_info_kappa = 5,
  factor_rank = 5,
  factor_features = 2000,
  factor_ridge = 1
)

fit$expression
fit$detection$membership_stats
fit$membership_fit$hierarchies
fit$local_geometry$W
fit$events[, c(
  "gene", "cell", "detection_block", "alra_margin", "recovered",
  "local_positive_mean", "local_positive_prevalence",
  "effective_donors", "predictability", "shrinkage",
  "prediction_sd"
)]
fit$predictive_variance
validate_dropout_result(fit, x)
```

## Direct recovery for a trusted mask

```r
membership_fit <- build_supercell_membership(
  embedding = pca,
  group = major_cell_type,
  split_by = condition,
  gamma = 150
)

rec <- recover_dropout_expression(
  x = x,
  mask = dropout_mask,
  membership = membership_fit,
  embedding = pca,
  recovery_method = "tree_local_factor",
  factor_target = "positive",
  return_details = TRUE
)
```

Both a bare membership vector and a full `DropoutKillerMembership` object use the final membership as the recovery borrowing block. A supplied hard stratum can further split that block but never expands borrowing beyond the final membership.

## Event diagnostics

Global-ALRA detection events include:

- `detection_block`
- `lowrank`
- `threshold`
- `alra_margin = lowrank - threshold`

Tree-local recovery adds:

- `local_positive_mean`
- `local_positive_variance`
- `local_positive_prevalence`
- `effective_donors` (Kish size of the **positive baseline** donor weights)
- `n_observed_gene` / `n_donors` (residual-model donor count available to the query)
- `tree_distance_weighted_mean`
- `embedding_distance_weighted_mean`
- `factor_prediction`
- `predictability`
- `shrinkage`
- `prediction_sd`
- `recovery_method` (`tree_local_mean`, `tree_local_factor`, or `unavailable`)

## Required recovery benchmark

Recovery must be validated with the same oracle masks across component ablations:

1. positive membership mean;
2. embedding-only kernel mean;
3. tree-only weighted mean;
4. tree + embedding weighted mean;
5. current `masked_factor`;
6. local mean + residual factor;
7. batched tree-local residual factor.

Required scenarios include MCAR positive masking, original-UMI count strata, Binomial-zero, and full binomial thinning followed by re-normalization, PCA, SuperCell reconstruction, detection, and recovery.

Report RMSE, MAE, bias, Pearson/CCC, interval coverage, gene-variance error, distributional distance, tree-distance strata, effective-donor strata, and biological-state strata. A post-recovery increase in correlation is not independent validation because coexpression participates in prediction.

## Comparator engines

Historical membership-local EB detector:

```r
fit_eb <- dropout_killer(
  x = x,
  embedding = pca,
  group = major_cell_type,
  detection_method = "eb_zero_null",
  threshold = 0.95
)
```

Historical membership-local ALRA-quantile detector:

```r
fit_local_alra <- dropout_killer(
  x = x,
  embedding = pca,
  group = major_cell_type,
  detection_method = "alra_quantile"
)
```

Current masked-factor recovery comparator:

```r
fit_masked_factor <- dropout_killer(
  x = x,
  embedding = pca,
  group = major_cell_type,
  recovery_method = "masked_factor",
  factor_target = "positive"
)
```

Historical Gaussian neighbor recovery comparator:

```r
fit_neighbor <- dropout_killer(
  x = x,
  embedding = pca,
  membership = membership,
  recovery_method = "neighbor"
)
```

The neighbor engine has no calibrated predictive-variance model, so `uncertainty_available = FALSE` and uncertainty-aware sampling rejects it.

## Seurat workflow

```r
obj <- dropout_killer_seurat(
  object = obj,
  assay = "RNA",
  slot = "counts",
  reduction = "pca",
  dims = 1:20,
  group_by = "major_cell_type",
  split_by = "condition",
  new_assay = "DropoutKiller"
)
```

With the default detector, `group_by` should identify the major cell class used for global ALRA zero detection. `split_by` remains a recovery boundary. Recovered values are continuous normalized-expression estimates and belong in assay data, not raw integer counts. Do not treat the completed mean matrix as an error-free count matrix for DE, trajectory, or network inference.

See `inst/ALGORITHM.md` for the historical detector/v0.5 factor contract and `inst/TREE_LOCAL_RECOVERY.md` for the v0.6 hierarchy-aware recovery contract.
