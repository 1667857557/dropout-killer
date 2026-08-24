# DropoutKiller

`DropoutKiller` is an R package for **selective** scRNA-seq dropout recovery. It never rewrites observed non-dropout coordinates and does not use external PPI/pathway/GRN priors.

Version 0.5 separates the problem explicitly into two statistical stages:

```text
expression X + cell embedding + optional hard biological strata
                         |
                         v
       SuperCell-style memberships (gamma = 150)
                         |
                         v
       membership-local low-rank reconstruction
                         |
                         v
 empirical-Bayes symmetric biological-zero null
 negative-tail variance shrinkage + gene-wise BH
                         |
                         v
             selective dropout mask M
                         |
                         v
    membership-local non-target coexpression factors
       (all current target genes are excluded)
                         |
                         v
 target-gene positive-conditional ridge prediction
                         |
                         v
       analytic leave-one-out optimal shrinkage
                         |
              +----------+----------+
              |                     |
              v                     v
       recovered mean        predictive variance
              |                     |
              +----------+----------+
                         v
       selective mean matrix + DV-aware draws
```

## Why the detector changed

The historical ALRA-derived detector used a per-gene empirical lower quantile with `quantile_prob = 0.001`. In a membership with `n` cells, that tail contains only `0.001 n` expected observations. For memberships below roughly 1000 cells, the empirical 0.1% quantile is therefore essentially determined by the first one or two order statistics and can be highly unstable.

The default detector now keeps the useful ALRA idea that biological-zero low-rank reconstruction error is approximately symmetric around zero, but estimates a finite-sample null instead of an extreme empirical quantile.

For gene `g` in one membership, let negative reconstructed values be `z_gc < 0`. Their second moment estimates the zero-null variance:

```text
s_g^2 = sum_{c:z_gc<0} z_gc^2 / n_g^-
```

Because `n_g^-` can be small, this is shrunk toward a robust membership-level variance center `s_0^2`:

```text
w_g = n_g^- / (n_g^- + nu_0)
s_g,EB^2 = w_g s_g^2 + (1-w_g) s_0^2
```

where `nu_0 = variance_prior_df` is the equivalent prior negative-residual sample size.

For an observed zero with positive low-rank reconstruction `z_gc`, the working biological-zero null gives

```text
Z_gc = z_gc / s_g,EB
p_gc = P(N(0,1) >= Z_gc)
```

P-values are Benjamini-Hochberg adjusted **within each gene and membership**. Event score is

```text
confidence_gc = 1 - q_gc
```

and is not a Bayesian posterior probability. With the default `threshold = 0.95`, selected events therefore satisfy approximately `q_gc <= 0.05` under the working symmetric zero-null model.

The historical empirical-quantile gate is still available with:

```r
detection_method = "alra_quantile"
```

for reproducibility.

## Recovery estimand

Once an event has already passed the technical-dropout detector, the recovery target changes. The desired magnitude is no longer the unconditional membership mean containing biological zeros. The default factor engine estimates

```text
E[X_gc | X_gc > 0, cell state, membership]
```

from reliable positive donor cells for target gene `g`.

Unmasked zeros remain in the original matrix and are never modified; they are simply not used to estimate the magnitude of a target that has already been classified as technical dropout. Set

```r
factor_target = "all_observed"
```

to reproduce the previous zero-inflated target model.

Inside a membership, every current recovery-target gene is excluded from the factor-feature set, so target expression cannot leak back into the cell-state representation used to predict it. Factor state is learned from standardized high-variance **non-target genes** by truncated SVD (`irlba` on large matrices).

For target gene `g`, reliable donor values are fit by ridge regression:

```text
beta_g = argmin_beta ||x_g,donor - X_donor beta||^2 + lambda ||beta_factor||^2
```

The ridge smoother provides an exact analytic leave-one-out prediction:

```text
yhat_i^(-i) = y_i - (y_i - yhat_i) / (1 - h_ii)
```

Let `mu_i^(-i)` be the leave-one-out positive-donor mean and

```text
d_i = yhat_i^(-i) - mu_i^(-i)
t_i = y_i - mu_i^(-i)
```

The factor contribution is shrunk by the held-out squared-error optimum:

```text
q_g = clip(sum_i d_i t_i / sum_i d_i^2, 0, 1)
Xhat_gc = max(0, mu_g + q_g (Xfactor_gc - mu_g))
```

Thus unsupported coexpression cannot force a cell-specific value: `q_g = 0` collapses to the positive membership mean.

Predictive variance is based on leave-one-out residual MSE rather than in-sample residual variance. For the default positive target, `sample_dropout_expression()` draws from a lower-truncated Gaussian approximation so a selected technical dropout is not sampled back into the biological-zero state.

## Matrix workflow

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
  recovery_method = "masked_factor",
  factor_target = "positive",
  factor_rank = 5,
  factor_features = 2000,
  factor_ridge = 1
)

fit$expression
fit$events[, c("gene", "cell", "q_value", "recovered",
               "prediction_sd", "predictability", "shrinkage")]
fit$predictive_variance
validate_dropout_result(fit, x)
```

If condition-specific states must never borrow from one another, pass condition to `split_by`. If cross-condition borrowing inside a trusted shared state is intended, leave `split_by = NULL`.

## Direct recovery for a supplied dropout mask

If the mask is already trusted, detection can be bypassed:

```r
rec <- recover_dropout_expression(
  x = x,
  mask = dropout_mask,
  membership = membership,
  recovery_method = "masked_factor",
  factor_target = "positive",
  return_details = TRUE
)
```

`embedding` is not required by the masked-factor engine.

Historical positional slots are retained. New 0.5 controls are appended after the 0.4 public API.

## Differential variability

Replacing missing values by conditional means alone contracts variance. For latent expression,

```text
Var(lambda_g | Y)
  = Var_i(E[lambda_ig | Y])
  + E_i(Var[lambda_ig | Y])
```

The deterministic mean matrix contains only the first component. The masked-factor engine therefore stores event-level `prediction_sd` and sparse `predictive_variance`. For DV, covariance, or coexpression analysis, propagate uncertainty with repeated completed draws:

```r
draws <- sample_dropout_expression(fit, n = 20, seed = 1)
```

Observed coordinates are identical in every draw.

The completed mean matrix should not be treated as an error-free raw-count matrix for DE, trajectory, or network inference. Recovered values are continuous conditional-expression estimates.

## Validation contract

A larger post-recovery correlation is not independent evidence that recovery was correct, because the same coexpression structure was used to predict the missing value. Validation should use held-out data:

1. pseudo-mask reliable observed positives and test predictive likelihood/error;
2. preferably perform count-level binomial/Poisson thinning so newly created zeros have known technical origin;
3. compare DV/covariance on repeated draws rather than only the mean matrix;
4. evaluate detection recall and mask stability separately from recovery accuracy.

## Legacy/comparator engines

Historical detection:

```r
fit_old_detect <- dropout_killer(
  x, pca,
  detection_method = "alra_quantile",
  quantile_prob = 0.001
)
```

Historical Gaussian neighbor recovery:

```r
fit_neighbor <- dropout_killer(
  x = x,
  embedding = pca,
  membership = membership,
  recovery_method = "neighbor",
  neighbor_k = 30
)
```

`weighted_neighbor_prediction()` is retained. Positive-only neighbor borrowing estimates a local conditional-positive mean, so it is a useful comparator when the masked coordinate is known to have been positive, but it must not be generalized to arbitrary natural zeros.

The neighbor engine has no calibrated predictive-variance model; `uncertainty_available` is therefore `FALSE` and `sample_dropout_expression()` rejects neighbor results.

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

Recovered values are continuous and are stored as assay data, not raw counts.

## Main API

- `build_supercell_membership()`: graph-based membership construction with canonical labels.
- `local_alra_detect()`: low-rank zero detection; finite-sample EB null by default.
- `local_alra_score()`: sparse `1 - q` evidence matrix.
- `select_dropout_mask()`: sparse high-confidence zero mask.
- `masked_factor_prediction()`: target-leakage-free membership-local coexpression prediction.
- `recover_dropout_expression()`: selective recovery for a supplied mask.
- `sample_dropout_expression()`: uncertainty-aware completed-matrix draws.
- `weighted_neighbor_prediction()`: legacy/comparator neighbor prediction.
- `dropout_killer()`: end-to-end workflow.
- `dropout_killer_seurat()`: Seurat wrapper.
- `validate_dropout_result()`: selective-recovery invariants.

See `inst/ALGORITHM.md` for the mathematical contract.
