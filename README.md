# DropoutKiller

`DropoutKiller` is an R package for **selective** scRNA-seq dropout recovery. It detects high-confidence zero events, then modifies only those coordinates. External PPI/pathway/GRN priors are not used.

## Workflow

```text
expression X + cell embedding + optional biological strata
                         |
                         v
       SuperCell-style memberships (gamma = 150)
                         |
                         v
        membership-local ALRA zero detection
                         |
                         v
             high-confidence dropout mask M
                         |
                         v
     masked membership-local coexpression factors
   (target genes excluded; M entries are missing)
                         |
                         v
       target-gene ridge regression on cell factors
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

## Recovery model

Inside membership `m`, let `R_gc = 0` only for dropout-mask coordinates and `R_gc = 1` otherwise. Every current recovery-target gene is excluded from the factor-feature set, so target expression cannot leak back into the cell-state representation used to predict it. The factor stage learns the remaining predictable coexpression component from non-target genes.

The implementation uses iterative masked low-rank reconstruction on standardized high-variance features to estimate cell factors `z_c`. For target gene `g`, only cells whose target value is not masked are used in ridge regression:

```text
beta_g = argmin_beta ||x_g,obs - X_obs beta||^2 + lambda ||beta_factor||^2
```

For fixed ridge penalty, the diagonal of the smoother matrix gives an exact analytic leave-one-out prediction:

```text
yhat_i^(-i) = y_i - (y_i - yhat_i) / (1 - h_ii)
```

Let `mu_i^(-i)` be the leave-one-out membership mean and

```text
d_i = yhat_i^(-i) - mu_i^(-i)
t_i = y_i - mu_i^(-i)
```

The shrinkage coefficient used for recovery is the squared-error optimum on these held-out predictions:

```text
q_g = clip(sum_i d_i t_i / sum_i d_i^2, 0, 1)
Xhat_gc = max(0, mu_g + q_g (Xfactor_gc - mu_g))
```

Thus unsupported coexpression cannot force a cell-specific value: `q_g = 0` is always available and collapses recovery to the membership mean. This requires no explicit CV folds. A separate leave-one-out predictability statistic is stored to summarize improvement over the membership-mean null.

Importantly, **unmasked zeros remain in the target-gene training data**; the model does not estimate `E[X | X>0]`.

A residual predictive variance is also retained. The deterministic recovered matrix contains predictive means, while `fit$predictive_variance` and `fit$events$prediction_sd` carry uncertainty needed for differential-variability-aware analysis.

## Matrix workflow

```r
library(DropoutKiller)

fit <- dropout_killer(
  x = x,
  embedding = pca,
  group = major_cell_type,
  split_by = condition,
  gamma = 150,
  recovery_method = "masked_factor",
  factor_rank = 5,
  factor_features = 2000,
  factor_ridge = 1
)

fit$expression
fit$events[, c("gene", "cell", "recovered", "prediction_sd", "predictability")]
fit$predictive_variance
validate_dropout_result(fit, x)
```

If condition-specific states must never borrow from one another, pass condition to `split_by`. If cross-condition borrowing inside a trusted shared state is intended, leave `split_by = NULL`.

## Direct recovery for a supplied dropout mask

`embedding` is not required by the default masked-factor engine:

```r
rec <- recover_dropout_expression(
  x = x,
  mask = dropout_mask,
  membership = membership,
  recovery_method = "masked_factor",
  return_details = TRUE
)
```

## Differential variability

Replacing missing values by conditional means alone necessarily contracts variance. For a latent expression value,

```text
Var(lambda_g | Y) = Var_i(E[lambda_ig | Y]) + E_i(Var[lambda_ig | Y])
```

The point matrix represents the first term. DropoutKiller therefore keeps the second term as event-level predictive variance. For downstream DV/coexpression analyses, use multiple completed draws rather than treating the mean matrix as error-free:

```r
draws <- sample_dropout_expression(fit, n = 20, seed = 1)
```

Observed entries are identical in every draw; only recovered dropout events vary.

## Legacy/comparator neighbor engine

The previous Gaussian neighbor estimator remains available explicitly:

```r
fit_neighbor <- dropout_killer(
  x = x,
  embedding = pca,
  membership = membership,
  recovery_method = "neighbor",
  neighbor_k = 30
)
```

`weighted_neighbor_prediction()` is also retained. Positive-only borrowing remains available for reproducibility, but its estimate is the local conditional-positive mean and can be upward-biased for dropout recovery.

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

- `build_supercell_membership()`: graph-based membership construction.
- `local_alra_detect()`: membership-local ALRA-inspired zero detection.
- `select_dropout_mask()`: sparse high-confidence zero mask.
- `masked_factor_prediction()`: masked membership-local coexpression prediction.
- `recover_dropout_expression()`: selective recovery for a supplied mask.
- `sample_dropout_expression()`: uncertainty-aware completed-matrix draws.
- `weighted_neighbor_prediction()`: legacy/comparator neighbor prediction.
- `dropout_killer()`: end-to-end workflow.
- `dropout_killer_seurat()`: Seurat wrapper.
- `validate_dropout_result()`: selective-recovery invariants.

See `inst/ALGORITHM.md` for the mathematical contract.

## Statistical boundaries

- The dropout detector and recovery model are separate. Recovery never changes the mask.
- The mean completed matrix is not claimed to preserve full DV by itself; use predictive variance or repeated draws for DV-sensitive downstream work.
- Target genes are excluded from factor-state learning to prevent target-to-predictor leakage.
- Factor rank controls only the predictable coexpression component. Residual variation is retained separately rather than forced into the low-rank mean.
- A gene with no held-out predictive gain shrinks toward its membership mean instead of receiving an unsupported coexpression estimate.
- Recovery is membership-local and never uses PPI/pathway/GRN priors.
