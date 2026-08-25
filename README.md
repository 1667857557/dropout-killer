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
                                     membership-local donor kernel
                                                       |
                                      query-specific positive mean
                                                       |
                                   gene-level coexpression residual
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

For production recovery the final membership is the computational and biological borrowing block:

- **explicit `group` / `split_by`**: additional absolute biological boundary;
- **final SuperCell membership**: absolute recovery borrowing block;
- **walktrap hierarchy + original embedding**: continuous weighting **inside** that membership.

This matches the intended use of `gamma`: first restrict recovery to a small biologically coherent cell set, then let closer cells contribute more than distant cells without repeatedly searching the whole hard stratum.

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

Thus millions of high-confidence dropout events increase output/indexing work but no longer create millions of independent regression fits.

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

Both a bare membership vector and a full `DropoutKillerMembership` object use the final membership as the recovery borrowing block. A supplied hard stratum can further split that block but never expands borrowing beyond the final membership.

## Event diagnostics

Tree-local events include:

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

## Required benchmark before default promotion

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