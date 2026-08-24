# DropoutKiller mathematical and engineering contract

## 1. Scope

Let `X in R_+^(G x C)` be the input expression matrix and

`D = {(g,c): X_gc = 0 and evidence supports technical dropout}`.

Only coordinates in `D` are eligible for replacement. Define

`R_gc = 0` for `(g,c) in D`, and `R_gc = 1` otherwise.

The key recovery contract is that `R_gc=0` means **missing**, not observed zero. Every unmasked zero remains an observed biological/sampling zero and participates in target-gene fitting.

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

while keeping unmasked zeros in the target-gene training data.

## 5. Leakage-free membership-local factor state

Within membership `k`, define the set of all genes currently targeted for recovery:

`T_k = {g : exists c with (g,c) in D_k}`.

Genes in `T_k` are excluded from factor-feature learning. This is deliberately conservative: the target gene and other simultaneously recovered target genes cannot contribute to the cell-state representation used to predict them.

From genes outside `T_k`, select up to `F` high-variance features with enough unmasked observations. For each selected gene, compute mean and variance from its observed entries and standardize

`Z_gc = (X_gc - mu_g) / s_g`.

Any masked coordinates among candidate features would be treated as missing, but because all current recovery targets are excluded, the default factor matrix is typically fully observed with respect to the active dropout mask. Missing working values are initialized at standardized mean zero. Iterative masked low-rank reconstruction updates only missing working coordinates:

`Z_mis^(t+1) = P_r(Z_obs + Z_mis^(t))_mis`,

where `P_r` is the rank-`r` orthogonal low-rank projection. Observed coordinates are never replaced during factor learning.

The resulting right singular/eigen vectors provide membership-local cell factors

`z_c in R^r`.

Because this shared factor stage uses an `n_k x n_k` Gram matrix after feature selection, it avoids constructing a dense `G x G` gene-correlation matrix.

## 6. Target-gene ridge prediction

For each target gene `g`, let `O_g` be membership cells whose `(g,c)` value is not masked. Unmasked zeros belong to `O_g` exactly like any other observed value.

Fit

`x_g,O = beta_0 + Z_O beta_g + epsilon`

by ridge regression:

`beta_hat = argmin_beta ||x_g,O - X_O beta||_2^2 + lambda ||beta_factor||_2^2`.

The intercept is unpenalized. If `P` is the ridge penalty matrix,

`beta_hat = (X_O^T X_O + P)^(-1) X_O^T x_g,O`.

The associated linear smoother is

`H = X_O (X_O^T X_O + P)^(-1) X_O^T`.

## 7. Exact analytic leave-one-out shrinkage

For a fixed ridge penalty, leave-one-out predictions can be obtained without explicitly refitting `n_g` models. Let `h_ii` be the diagonal of `H`, `yhat_i` the full-data fitted value, and `e_i = y_i-yhat_i`. Then

`yhat_i^(-i) = y_i - e_i/(1-h_ii)`.

Let the leave-one-out membership-mean null be

`mu_i^(-i) = (sum_j y_j - y_i)/(n_g-1)`.

Define the held-out factor deviation and held-out target deviation:

`d_i = yhat_i^(-i) - mu_i^(-i)`,

`t_i = y_i - mu_i^(-i)`.

Consider the convex combination

`ytilde_i(q) = mu_i^(-i) + q d_i`, `0 <= q <= 1`.

Its leave-one-out squared error is

`L(q) = sum_i (t_i - q d_i)^2`.

The unconstrained minimizer follows directly from `dL/dq=0`:

`q_raw = sum_i d_i t_i / sum_i d_i^2`.

The recovery shrinkage is therefore

`q_g = clip(q_raw, 0, 1)`.

This is the exact squared-error optimum between the ridge factor predictor and membership-mean null on analytic leave-one-out predictions. Because `q=0` is always admissible, unsupported coexpression cannot force a cell-specific prediction.

For a masked query cell `c`, let `m_gc^factor` be its full ridge factor prediction. The deterministic recovered mean is

`m_gc = max[0, mu_g + q_g (m_gc^factor - mu_g)]`.

A separate held-out predictability statistic is

`P_g = clip(1 - SSE_LOO,shrunk/SSE_LOO,null, 0, 1)`.

`q_g` controls shrinkage; `P_g` summarizes observed leave-one-out predictive gain. Small or under-supported target fits fall back to `q_g=0` and the membership mean.

## 8. Predictive uncertainty and differential variability

A conditional mean alone cannot preserve full differential variability. Let

`df_g = tr[(X_O^T X_O + P)^(-1)X_O^T X_O]`

and use the shrunken effective degrees of freedom

`df_g^* = 1 + q_g(df_g - 1)`.

Estimate residual variance from the shrunken in-sample mean:

`sigma_g^2 = sum_{c in O_g}(x_gc - m_gc^fit)^2 / max(n_g - df_g^*, 1)`.

For query design row `x_c`, ridge leverage is

`h_gc = x_c^T (X_O^T X_O + P)^(-1) x_c`.

The event-level predictive variance approximation is

`v_gc = sigma_g^2 [1 + q_g^2 h_gc]`.

DropoutKiller therefore represents each recovered event by both

`m_gc ~= E[X_gc | observed information]`

and

`v_gc ~= Var[X_gc | observed information]`.

For a gene, the variance decomposition motivating DV-aware downstream analysis is

`Var(X_g | Y) = Var_c(E[X_gc | Y]) + E_c(Var[X_gc | Y])`.

A point-imputed matrix retains only the first term. `predictive_variance` / `prediction_sd` retain the second term explicitly.

## 9. Uncertainty-aware completed draws

`sample_dropout_expression()` generates repeated completed matrices. For a recovered event,

`X_gc^(b) = max(0, Normal(m_gc, v_gc))`.

Observed coordinates are fixed exactly in every draw. The Gaussian predictive distribution is an approximation on the same expression scale as the package input; these draws are intended for sensitivity, DV, and coexpression propagation rather than reinterpretation as raw UMI counts.

Repeated draws allow downstream summaries such as

`E_b[Var_c(X_g^(b))]`

or repeated covariance/network estimation without treating every recovered mean as perfectly known.

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

For membership size `n_k`, selected non-target factor features `F`, factor rank `r`, and a small number of masked-factor iterations:

- factor state uses feature-by-cell operations plus an `n_k x n_k` Gram eigendecomposition;
- each target gene solves an `(r+1) x (r+1)` ridge system;
- exact LOO diagnostics require only the smoother diagonal, not `n_g` separate refits;
- no dense `G x G` coexpression matrix is formed.

With default `gamma=150`, `r=5`, and `F<=2000`, this is substantially smaller than a full graphical model or all-gene elastic-net inside every membership.

## 13. Invariants

1. Mask entries must correspond to original zeros.
2. Masked zeros are missing for recovery; unmasked zeros remain target observations.
3. Current recovery-target genes do not contribute to factor-state learning.
4. Observed non-dropout values remain numerically exact.
5. Recovery never crosses membership boundaries.
6. Unsupported factor predictions shrink to the membership mean.
7. Predictive uncertainty is retained separately from the deterministic mean matrix.
8. External biological priors do not enter detection or recovery.
9. The dropout mask is never redefined by recovery.
10. Iterative factor updates change only missing working coordinates, not observed data.

## 14. Validation principle

Recovery correctness cannot be established by observing stronger correlation after imputation, because the same coexpression structure generated the prediction. The analytic leave-one-out target regression provides an internal, leakage-reduced prediction check at negligible extra fitting cost. Larger empirical benchmarks should still use pseudo-dropout or Poisson/binomial thinning and evaluate held-out predictive error, DV recovery, and covariance recovery against independent information.
