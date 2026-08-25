# DropoutKiller

`DropoutKiller` is an R package for **selective** scRNA-seq dropout recovery. It never overwrites observed non-dropout coordinates and does not use external PPI/pathway/GRN priors.

Version 0.6 keeps detection and recovery as separate statistical problems and adds a hierarchy-aware recovery engine without replacing the validated v0.5 comparator by default.

```text
expression X + biological embedding + optional hard strata
                         |
                         v
          SuperCell kNN graph + walktrap hierarchy
                         |
                         +-------------------------+
                         |                         |
                         v                         v
                 final gamma cut            retained tree
                         |                         |
                         +------------+------------+
                                      |
                                      v
                        EB zero-null dropout detector
                                      |
                                      v
                             selective dropout mask
                                      |
                      +---------------+----------------+
                      |                                |
                      v                                v
             masked_factor                    tree_local_factor
              comparator                           v0.6 engine
                                                       |
                                      tree + embedding donor kernel
                                                       |
                                      query-specific positive mean
                                                       |
                                      local coexpression residual
                                                       |
                                      held-out shrinkage + variance
                                                       |
                                                       v
                                      selective means + Gamma draws
```

## Detection

The default detector remains `detection_method = "eb_zero_null"`. It replaces the old membership-local empirical `quantile_prob = 0.001` threshold with a finite-sample zero-null model.

For negative low-rank reconstructions of gene `g`,

```text
s_g^2 = sum(z_gc^2 : z_gc < 0) / n_g^-
w_g   = n_g^- / (n_g^- + variance_prior_df)
s_g,EB^2 = w_g s_g^2 + (1-w_g) s_0^2
```

and a positive reconstruction at an observed zero is tested by

```text
Z_gc = z_gc / s_g,EB
p_gc = P(N(0,1) >= Z_gc)
```

with gene-wise BH adjustment inside each membership. `confidence = 1 - q_value`; it is **not** a Bayesian posterior probability. The historical detector remains available as `detection_method = "alra_quantile"`.

## Why retain the SuperCell hierarchy?

A SuperCell membership is only a cut through a walktrap hierarchy:

```text
biological embedding -> kNN graph -> walktrap tree -> gamma cut
```

Treating all cells in the final cut as exchangeable discards information already present in that tree. Version 0.6 therefore keeps the `cluster_walktrap()` hierarchy and a compact tree index inside `DropoutKillerMembership`.

Hard biological strata and the final membership have different roles:

- **explicit `group` / `split_by`**: absolute borrowing boundary;
- **final SuperCell membership**: high-weight local core inside an explicit hard stratum, not an artificial zero-weight wall;
- **if no hard stratum is explicitly supplied**: the final membership remains the conservative borrowing boundary.

This prevents accidental cross-lineage borrowing in heterogeneous datasets while allowing nearby sibling memberships to contribute when the user has supplied a trusted biological stratum.

## Tree-local donor weighting

For query cell `c` and candidate donor `j` in the same hard stratum, the new engine uses

```text
w_cj = exp[-alpha d_tree(c,j)/tau
           -(1-alpha) d_embed(c,j)^2/(2 h_c^2)]
```

where:

- `alpha = tree_weight`;
- `d_tree` is a normalized walktrap lowest-common-ancestor distance;
- `d_embed` is distance in the original PCA/Harmony/WNN-like biological embedding;
- `h_c` is an adaptive bandwidth defined by the `local_k`-th candidate neighbor.

If hierarchy distance is unavailable for a pair (for example a non-anchor cell in an approximate build), the pair uses the embedding component instead of inventing a tree distance.

The sparse candidate graph is built once per recovery run and reused across genes.

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

## Local residual coexpression

The tree-local engine does not ask coexpression factors to predict the entire target expression. It first defines leave-one-out local residuals

```text
r_gj = X_gj - mu_gj^(-j)
```

and predicts only what the tree/embedding-local positive baseline cannot explain.

All current target genes are excluded from factor-state learning, preserving the no-target-leakage path

```text
X_cell,-targets -> factor state -> target residual
```

For query `c`, residual ridge is locally weighted:

```text
beta_gc = (F' W_c F + P)^-1 F' W_c r_g
```

so biologically nearer donors have greater influence on both the baseline and the coexpression correction.

Held-out residual predictions estimate `q_pred`. Effective residual-donor information contributes

```text
q_info = n_eff / (n_eff + local_info_kappa)
q_final = q_pred * q_info
```

and the deployed prediction is

```text
Xhat_gc = max(0, mu_gc + q_final * rhat_gc)
```

If the local factor model has no held-out gain, recovery falls back to the **query-specific local positive mean**, not the whole-membership mean.

## Predictive uncertainty and DV

Tree-local fallback variance retains both positive-expression outcome variance and local-mean estimation uncertainty. For a supported residual model, held-out error is recomputed using the same deployed `q_final` used in the recovered mean, rather than the stronger pre-information-shrinkage coefficient.

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

The v0.6 engine is introduced as a **parallel benchmark engine** first; `masked_factor` remains the default until independent oracle/full-thinning benchmarks establish superiority.

```r
library(DropoutKiller)

fit <- dropout_killer(
  x = x,
  embedding = pca,
  group = major_cell_type,
  split_by = condition,
  gamma = 150,
  detection_method = "eb_zero_null",
  variance_prior_df = 10,
  threshold = 0.95,
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
fit$membership_fit$hierarchies
fit$local_geometry$W
fit$events[, c(
  "gene", "cell", "q_value", "recovered",
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

If only a bare membership vector is supplied and no explicit hard stratum is provided, tree-local recovery conservatively stays within that final membership. To permit borrowing across sibling memberships, provide a trusted hard stratum via the full `DropoutKillerMembership` object or `hard_stratum` in direct recovery.

## Event diagnostics

Tree-local events include:

- `local_positive_mean`
- `local_positive_variance`
- `local_positive_prevalence`
- `effective_donors` (Kish size of the **positive baseline** donor weights)
- `n_observed_gene` / `n_donors` (residual-model donor count)
- `tree_distance_weighted_mean`
- `embedding_distance_weighted_mean`
- `factor_prediction`
- `predictability`
- `shrinkage`
- `prediction_sd`
- `recovery_method` (`tree_local_mean`, `tree_local_factor`, or `unavailable`)

## Required benchmark before default promotion

Recovery must be validated with the same oracle masks across component ablations:

1. positive membership mean;
2. embedding-only kernel mean;
3. tree-only weighted mean;
4. tree + embedding weighted mean;
5. current `masked_factor`;
6. local mean + residual factor;
7. full tree-local weighted residual factor.

Required scenarios include MCAR positive masking, original-UMI count strata, Binomial-zero, and full binomial thinning followed by re-normalization, PCA, SuperCell reconstruction, detection, and recovery.

Report RMSE, MAE, bias, Pearson/CCC, interval coverage, gene-variance error, distributional distance, tree-distance strata, effective-donor strata, and biological-state strata. A post-recovery increase in correlation is not independent validation because coexpression participates in prediction.

## Comparator engines

Current masked-factor comparator:

```r
fit_masked_factor <- dropout_killer(
  x = x,
  embedding = pca,
  group = major_cell_type,
  recovery_method = "masked_factor",
  factor_target = "positive"
)
```

Historical Gaussian neighbor comparator:

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
  slot = "data",
  reduction = "pca",
  dims = 1:20,
  group_by = "major_cell_type",
  split_by = "condition",
  new_assay = "DropoutKiller"
)
```

Recovered values are continuous normalized-expression estimates and belong in assay data, not raw integer counts. Do not treat the completed mean matrix as an error-free count matrix for DE, trajectory, or network inference.

See `inst/ALGORITHM.md` for the detector/v0.5 factor contract and `inst/TREE_LOCAL_RECOVERY.md` for the v0.6 hierarchy-aware recovery contract.
