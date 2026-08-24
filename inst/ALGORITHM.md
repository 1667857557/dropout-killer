# DropoutKiller mathematical and engineering contract

## 1. Scope

Let `X in R_+^(G x C)` be the input expression matrix and

`D = {(g,c): X_gc = 0 and evidence supports technical dropout}`.

Only coordinates in `D` are eligible for replacement. Define

`R_gc = 0` for `(g,c) in D`, and `R_gc = 1` otherwise.

The key recovery contract is that `R_gc=0` means **missing**, not observed zero. Every unmasked zero remains an observed biological/sampling zero and participates in model fitting.

No external PPI, pathway, GRN, or TF-target prior is used.

## 2. Membership construction

For a low-dimensional cell representation `z_c`, memberships are constructed separately inside supplied hard biological strata. Within stratum `s`, build a Euclidean kNN graph and apply walktrap clustering. With `n_s` cells and graining level `gamma`, the target is

`K_s = max(1, round(n_s/gamma))`.

The default `gamma=150` therefore targets memberships of order 150 cells while preserving hard biological boundaries.

## 3. Dropout detection

Detection remains independent of recovery. For membership `k`, compute an uncentered low-rank approximation

`Y_k ~= U_k Sigma_k V_k^T = Yhat_k`.

For gene `g`, an observed zero is eligible only if its low-rank value passes the ALRA-derived negative-tail gate. The negative reconstruction tail supplies a working null scale and a one-sided confidence score. The default final mask is

`M_gc = I(X_gc=0) I(ALRA_gate_gc=1) I(C_gc>=0.95)`.

`C_gc` is evidence against the working biological-zero reconstruction null, not a calibrated posterior dropout probability.

## 4. Why recovery is no longer positive-only neighbor averaging

The previous default estimated

`E[X_g | X_g > 0, local state]`

because zero-valued donors were excluded. In general,

`E[X_g | X_g > 0, state] >= E[X_g | state]`,

so using the conditional-positive mean for a known missing expression value introduces upward bias whenever the target distribution has nonzero zero mass.

The new default instead estimates

`E[X_gc | X_c,-g observed, same membership, M_gc=1]`

while keeping unmasked zeros in the training data.

## 5. Masked membership-local factor state

Within membership `k`, select up to `F` high-variance genes with enough unmasked observations. For each selected gene, compute the mean and variance using only `R_gc=1` entries and standardize

`Z_gc = (X_gc - mu_g) / s_g`.

Masked entries are initialized at standardized mean zero. Iterative masked low-rank reconstruction updates only the missing coordinates:

`Z_mis^(t+1) = P_r(Z_obs + Z_mis^(t))_mis`,

where `P_r` is the rank-`r` orthogonal low-rank projection. Observed coordinates are never replaced during this factor-learning step.

The resulting right singular/eigen vectors provide membership-local cell factors

`z_c in R^r`.

Because the factor stage is shared across genes, its dominant cost is based on a `n_k x n_k` Gram matrix after feature selection rather than a `G x G` gene-correlation matrix.

## 6. Target-gene conditional prediction

For each target gene `g`, let `O_g` be membership cells whose `(g,c)` value is not masked. Fit ridge regression

`x_g,O = beta_0 + Z_O beta_g + epsilon`,

with

`beta_hat = argmin_beta ||x_g,O - X_O beta||_2^2 + lambda ||beta_factor||_2^2`.

The intercept is not penalized. If `X_O` is the design matrix and `P` the ridge penalty matrix,

`beta_hat = (X_O^T X_O + P)^(-1) X_O^T x_g,O`.

The effective degrees of freedom are

`df_g = tr[(X_O^T X_O + P)^(-1) X_O^T X_O]`.

## 7. GCV predictability shrinkage

A flexible factor model should not be allowed to create cell-specific expression merely because a low-rank representation exists. Compare its analytic generalized cross-validation error with the intercept-only membership model:

`GCV_factor = (SSE_factor/n) / (1 - df_g/n)^2`,

`GCV_null = (SSE_null/n) / (1 - 1/n)^2`.

Define

`q_g = clip(1 - GCV_factor/GCV_null, 0, 1)`.

The recovered conditional mean is

`m_gc = max[0, mu_g + q_g (m_gc^factor - mu_g)]`.

Therefore:

- if factor structure has reproducible predictive value, `q_g > 0` and cell state contributes;
- if it does not outperform the membership mean, `q_g = 0` and the model collapses to `mu_g`;
- small or under-supported target fits also fall back to `mu_g` rather than extrapolating unstable coefficients.

This is an analytic shrinkage gate and does not require explicit cross-validation folds.

## 8. Predictive uncertainty and differential variability

The deterministic completed matrix stores `m_gc`, but a conditional mean alone cannot preserve full differential variability. After GCV shrinkage, estimate residual variance

`sigma_g^2 = sum_{c in O_g}(x_gc - m_gc^fit)^2 / max(n_g - df_g^*, 1)`,

where

`df_g^* = 1 + q_g(df_g - 1)`.

For query design row `x_c`, ridge leverage is

`h_gc = x_c^T (X_O^T X_O + P)^(-1) x_c`.

The event-level predictive variance approximation is

`v_gc = sigma_g^2 [1 + q_g^2 h_gc]`.

DropoutKiller returns both

`m_gc = E[X_gc | observed information]`

and

`v_gc ~= Var[X_gc | observed information]`.

For a gene, the variance decomposition motivating DV-aware downstream analysis is

`Var(X_g | Y) = Var_c(E[X_gc | Y]) + E_c(Var[X_gc | Y])`.

A point-imputed matrix retains only the first term. `predictive_variance` / `prediction_sd` retain the second term explicitly.

## 9. Uncertainty-aware completed draws

`sample_dropout_expression()` generates repeated completed matrices. For a recovered event,

`X_gc^(b) = max(0, Normal(m_gc, v_gc))`.

Observed coordinates are fixed exactly in every draw. The Gaussian predictive distribution is an approximation on the same expression scale as the package input; the function is intended for sensitivity, DV, and coexpression propagation rather than reinterpretation as raw UMI counts.

Repeated draws allow downstream summaries such as

`E_b[Var_c(X_g^(b))]`

or repeated covariance/network estimation without falsely treating every recovered mean as perfectly known.

## 10. Legacy neighbor engine

The previous Gaussian cell-neighbor estimator remains available with `recovery_method="neighbor"` and through `weighted_neighbor_prediction()`. It is no longer the default.

For latent distance `d_cj`, neighbor weights are

`w_cj proportional to exp(-d_cj^2/sigma_c^2)`.

Positive-only donor renormalization remains available for reproducibility, but it targets a conditional-positive mean and is therefore not the default selective-recovery estimator.

## 11. Selective replacement

Final deterministic expression is

`Xfinal_gc = X_gc` if `M_gc=0`,

and

`Xfinal_gc = m_gc` if `M_gc=1` and a finite positive prediction is available.

Observed non-dropout values are never overwritten.

## 12. Computational scaling

For membership size `n_k`, selected factor features `F`, factor rank `r`, and a small number of masked-factor iterations:

- factor state uses feature-by-cell operations plus an `n_k x n_k` Gram eigendecomposition;
- each target gene solves an `(r+1) x (r+1)` ridge system;
- no dense `G x G` coexpression matrix is formed.

With the package default `gamma=150`, `r=5`, and `F<=2000`, this is substantially smaller than fitting a full graphical model or all-gene elastic-net inside every membership.

## 13. Invariants

1. Mask entries must correspond to original zeros.
2. Masked zeros are missing for recovery learning; unmasked zeros remain observations.
3. Observed non-dropout values remain numerically exact.
4. Recovery never crosses membership boundaries.
5. Unsupported factor predictions shrink to the membership mean.
6. Predictive uncertainty is retained separately from the deterministic mean matrix.
7. External biological priors do not enter detection or recovery.
8. The dropout mask is never redefined by recovery.
9. Iterative factor updates change only masked working values, not observed data.

## 14. Validation principle

Recovery correctness cannot be established by observing that correlation becomes stronger after imputation, because the same coexpression structure generated the prediction. Valid benchmarking should hide known values, refit without them, and compare held-out predictive likelihood/error against the membership-mean baseline. The internal GCV shrinkage serves as a low-cost gene-wise prediction gate; larger empirical benchmarks should still use pseudo-dropout or thinning-based held-out validation.
