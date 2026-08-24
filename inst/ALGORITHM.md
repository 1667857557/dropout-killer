# DropoutKiller mathematical and engineering contract

## 1. Scope

Let `X in R_+^(G x C)` be a normalized scRNA-seq expression matrix. DropoutKiller estimates only a sparse event set

`D = {(g,c): X_gc = 0 and evidence supports technical dropout}`

and only those coordinates are eligible for replacement. Recovery uses expression and latent cell geometry only; no external biological prior enters detection or recovery.

## 2. SuperCell-style membership

For a low-dimensional cell representation `z_c`, memberships are constructed separately inside hard strata such as major cell type and, optionally, condition or donor.

Within stratum `s`, construct a Euclidean kNN graph `G_s=(V_s,E_s)` and apply walktrap hierarchical community detection. With `n_s` cells and graining level `gamma`, the requested number of memberships is

`K_s = max(1, round(n_s/gamma))`.

The dendrogram is cut at `K_s`, or at least the number of disconnected graph components. For large strata an anchor subset is clustered first; omitted cells are assigned to the nearest membership centroid in embedding space.

## 3. Membership-local low-rank model

For membership `k`, let

`Y_k = X[, C_k]`.

Compute an uncentered rank-`r_k` approximation

`Y_k ~= U_k Sigma_k V_k^T = Yhat_k`.

When genes greatly outnumber cells, the implementation may use the equivalent Gram decomposition

`Y_k^T Y_k = V Lambda V^T`

and reconstruct

`Yhat_k = Y_k V_r V_r^T`.

A numeric rank can be supplied. With `rank="auto"`, singular-value spacings are compared with their tail distribution following the ALRA rank-selection principle.

## 4. ALRA-derived zero-preserving gate

For gene `g` in membership `k`, define

`q_gk = Q_p(Yhat_g,Ck)`

and, when the lower quantile is negative,

`tau_gk = |q_gk|`.

A value is an eligible candidate only if

`X_gc = 0` and `Yhat_gc > tau_gk`.

Observed non-zero expression is therefore excluded before recovery.

## 5. Confidence against the biological-zero null

Negative low-rank values estimate a symmetric reconstruction-error distribution associated with biological zeros. The working null is

`E_gk ~ N(0, sigma_gk^2)`

with

`sigma_gk^2 = mean(e^2 : e = Yhat_gj < 0)`.

For an ALRA-gated zero,

`C_gc = Phi(Yhat_gc / sigma_gk)`.

If too few negative values are available to estimate `sigma_gk`, confidence is set to `0.5`. The default mask is

`M_gc = I(X_gc=0) I(Yhat_gc>tau_gk) I(C_gc>=0.95)`.

`C_gc` is a null-tail confidence score, not a calibrated posterior dropout probability.

## 6. Estimated neighbor weights

For a masked event `(g,c)`, only cells in the same membership are eligible donors. With latent distance

`d_cj^2 = ||z_c-z_j||_2^2`,

raw Gaussian weights are

`a_cj = exp(-d_cj^2 / sigma_c^2)`.

If `neighbor_sigma=NULL`, the bandwidth is estimated for each query cell as the median positive distance among selected neighbors:

`sigma_c = median({d_cj : d_cj > 0})`.

The all-neighbor normalized weights are

`w_cj = a_cj / sum_l a_cl`.

By default, donors with `X_gj=0` are excluded for the target gene and the effective weights are renormalized over positive donors:

`w*_gcj = w_cj I(X_gj>0) / sum_l w_cl I(X_gl>0)`.

The recovered value is therefore

`Xhat_gc = sum_j w*_gcj X_gj`.

Setting `neighbor_positive_only=FALSE` instead uses

`Xhat_gc = sum_j w_cj X_gj`.

If the required donor support is unavailable, the zero remains unchanged. An optional `cap_quantile` can cap the prediction at a membership-local positive-expression quantile for the target gene.

## 7. Selective replacement

The final matrix is

`Xfinal_gc = X_gc` for `M_gc = 0`,

and

`Xfinal_gc = Xhat_gc` for `M_gc = 1` when a finite positive neighbor prediction is available.

There is no PPI/pathway branch and no fixed hybrid mixing coefficient.

## 8. Invariants

1. Mask entries must correspond to original zeros.
2. Observed non-zero values are unchanged exactly.
3. Recovery never crosses membership boundaries.
4. Supplied hard biological strata cannot be crossed during membership graph construction.
5. Recovery depends only on membership-local cell expression and estimated latent-space neighbor weights.
6. Recovered values are continuous expression estimates, not raw counts.
7. Sparse input remains sparse in the end-to-end recovery path.
8. Recovery is one-pass; recovered values are never fed back into detection or membership construction.

## 9. Statistical trade-offs

- Hard strata prevent cross-label borrowing but make recovery dependent on annotation quality.
- The Gaussian negative-tail model is a working null; its confidence is not an FDR or posterior probability.
- Positive-only borrowing targets the conditional-positive local mean and may be upward-biased.
- Very small memberships or weak positive donor support can leave eligible zeros unrecovered rather than forcing an estimate.
- The pipeline deliberately avoids iterative re-imputation because iterative reuse can amplify the package's own predictions.
