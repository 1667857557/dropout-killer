# DropoutKiller

`DropoutKiller` is an R package for **selective** scRNA-seq dropout recovery. Observed coordinates are never overwritten. Detection and recovery are separate statistical problems.

Version 0.8 changes only the production **dropout detector**. The default detector is now a source-faithful R implementation of the Gamma-Normal detector released with **scGACL**. The validated recovery default remains **P1_STABILIZED_STATE**.

```text
raw counts X + major cell class + biological embedding
        |
        +---------------- detector ----------------+
        |                                          |
        v                                          |
  log(1.01 + X)                                   |
        |                                          |
        v                                          |
 gene-wise Gamma-Normal EM                         |
 inside each major cell class                      |
        |                                          |
        v                                          |
 Gamma-component posterior >= 0.5                  |
        |                                          |
        +-------------- dropout mask --------------+
                                                   |
raw counts X                                       |
        |                                          |
        v                                          |
ALRA library-size + log1p recovery working scale   |
        |                                          |
        v                                          |
 SuperCell hierarchy + embedding                   |
        |                                          |
        v                                          |
 P1 stabilized-state selective recovery <----------+
```

## Default detection: scGACL Gamma-Normal mixture

The high-level default is:

```r
detection_method = "scgacl_gamma_normal"
scgacl_dropout_threshold = 0.5
```

This implementation follows the released scGACL detector code in [`Hhjyl/scGACL/evaluation/dropout_identify.py`](https://github.com/Hhjyl/scGACL/blob/main/evaluation/dropout_identify.py) and its released configuration (`IDENTIFY_USE_CELLTYPE = True`, `DROP_THRESHOLD = 0.5`).

### Detector input scale

The detector receives the **raw count matrix**, not the ALRA-normalized recovery matrix. For every raw count:

```text
x = log(1.01 + raw_count)
point = log(1.01)
```

Therefore every observed raw zero maps exactly to `point`.

`group` supplies the scGACL subpopulation labels. The released scGACL configuration uses known cell-type labels, so the source-faithful DropoutKiller path requires `group`. It does not silently substitute a different R clustering implementation when labels are absent.

### Per-gene mixture

For every gene inside every supplied group, scGACL fits

```text
f(x) = lambda Gamma(x; alpha, beta)
       + (1-lambda) Normal(x; mu, sigma)
```

where `beta` is the Gamma **rate**. The posterior responsibility of the Gamma component is

```text
d(x) = lambda Gamma(x; alpha, beta)
       / [lambda Gamma(x; alpha, beta)
          + (1-lambda) Normal(x; mu, sigma)]
```

The released initialization and EM contract are retained:

- `lambda` starts at the fraction of entries equal to `point`;
- genes with initial `lambda > 0.95` are invalid;
- if `lambda == 0`, it is initialized to `0.01`;
- `alpha = 1.5`, `beta = 1` initially;
- `mu` and population-SD `sigma` are initialized from values `x > point`;
- the Gamma update solves `log(alpha) - digamma(alpha) = v` and then sets `beta = alpha * sum(w) / sum(w*x)`;
- EM stops when the squared change in base-10 log likelihood is `<= 0.5`, or after the released iteration limit;
- genes whose transformed mean lies within `1e-2` of `point` are marked invalid before EM.

Only observed raw zeros are eligible for the DropoutKiller mask. The released scGACL code keeps a zero when its dropout posterior is not below the configured threshold; DropoutKiller therefore uses `d >= 0.5` at the exact boundary.

For a fixed gene and group, all raw zeros have the same transformed value `point`; evaluating the zero posterior once and assigning it to all zero coordinates is algebraically identical to constructing the full dense posterior matrix.

## Detection and recovery use different scales

This separation is intentional:

```text
scGACL detection: raw counts -> log(1.01 + count)
P1 recovery:       raw counts -> library normalize to 1e4 -> log1p
```

Changing `normalization_scale_factor` changes the recovery working scale but does **not** change the scGACL dropout mask.

## Optional original ALRA detector

ALRA remains available explicitly:

```r
detection_method = "alra_global"
```

This comparator follows the released KlugerLab/ALRA implementation on **one all-cell matrix**. `group`, `split_by`, and final SuperCell memberships do not fragment ALRA detection.

With normalized cell-by-gene matrix `A`, automatic rank selection uses up to the original `K = 100` singular values:

```text
d_i = sigma_i - sigma_{i+1}
```

The noise spacing distribution uses the original `noise_start = 80` tail and rank is the largest spacing more than six standard deviations above that tail mean. The original source validity checks are retained rather than adapting `K` or the noise tail to small cell-class blocks.

The final randomized rank-`k` SVD uses the original `q = 10`, followed by the original gene-wise zero threshold:

```text
tau_g = |Q_0.001(Ahat_.g)|
call_gc = (A_gc == 0) & (Ahat_gc > tau_g)
```

There is no second confidence gate for ALRA.

The old name:

```r
detection_method = "alra_global_by_group"
```

is deprecated. It emits a warning and delegates to `alra_global`; it no longer fits separate cell-class ALRA models.

Historical membership-local comparison engines remain available:

```r
detection_method = "eb_zero_null"
detection_method = "alra_quantile"
```

## Production recovery: P1_STABILIZED_STATE

The recovery default remains:

```r
recovery_method = "p1_stabilized_state"
```

For every deterministic target-gene fold, the target fold is excluded before predictor construction. Standardized non-target predictors receive one row-stochastic hierarchy/embedding smoothing step:

```text
Z_stable = (1-rho) Z + rho Z P'
```

with `rho = 0.25` by default. P1 then fits positive-donor ridge models with target-safe factors, analytic leave-one-out shrinkage, support-adaptive rank, and bias calibration. Final expression values are never graph-smoothed.

The SuperCell final membership remains an absolute recovery borrowing block. The retained walktrap hierarchy and original biological embedding provide continuous donor weighting inside that membership.

## Matrix workflow

```r
library(DropoutKiller)

fit <- dropout_killer(
  x = raw_counts,
  embedding = pca,
  group = major_cell_type,
  split_by = condition,
  gamma = 150,
  detection_method = "scgacl_gamma_normal",
  scgacl_dropout_threshold = 0.5,
  recovery_method = "p1_stabilized_state",
  factor_target = "positive",
  factor_rank = 5,
  factor_features = 2000,
  factor_ridge = 2,
  min_target_observed = 8,
  factor_crossfit_folds = 5,
  support_adaptive_rank = TRUE,
  bias_kappa = 10,
  predictor_smoothing = 0.25
)

fit$expression
fit$mask
fit$detection$membership_stats
fit$events[, c(
  "gene", "cell", "detection_block", "confidence", "recovered",
  "factor_rank", "prediction_sd", "recovery_method"
)]
```

For the original ALRA comparator:

```r
fit_alra <- dropout_killer(
  x = raw_counts,
  embedding = pca,
  membership = membership,
  detection_method = "alra_global",
  recovery_method = "p1_stabilized_state"
)
```

## scGACL event diagnostics

Selected scGACL events include:

- `detection_block`: supplied cell subpopulation;
- `confidence`: Gamma-component dropout posterior;
- `mixture_rate`;
- `gamma_shape`;
- `gamma_rate`;
- `normal_mean`;
- `normal_sd`.

Original-ALRA events instead retain:

- `lowrank`;
- `threshold`;
- `alra_margin = lowrank - threshold`.

P1 recovery adds recovery, donor, factor, shrinkage, and uncertainty diagnostics to either detector's event table.

## Direct recovery for a trusted mask

```r
rec <- recover_dropout_expression(
  x = normalized_expression,
  mask = trusted_dropout_mask,
  membership = membership_fit,
  embedding = pca,
  recovery_method = "p1_stabilized_state",
  factor_target = "positive",
  return_details = TRUE
)
```

## Seurat workflow

Use the raw `counts` slot and provide a major-cell-type metadata field for the default scGACL detector:

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

A pre-normalized Seurat `data` slot is not a source-faithful input for the default scGACL detector. Use raw counts, or explicitly select another detector whose input contract matches the supplied matrix.

## Method provenance

- scGACL: Huang et al., *Briefings in Bioinformatics* (2026), released code: <https://github.com/Hhjyl/scGACL>.
- ALRA: Linderman, Zhao & Kluger, *Nature Communications* 13, 192 (2022), released code: <https://github.com/KlugerLab/ALRA>.

DropoutKiller ports the detector mathematics; it does not vendor scGACL's GAN/VAE imputation model. P1 remains the DropoutKiller recovery engine after a detector has selected technical-zero coordinates.
