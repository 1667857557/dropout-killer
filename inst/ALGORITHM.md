# DropoutKiller mathematical and engineering contract

## 1. Scope

Let `X in R_+^(G x C)` be a normalized scRNA-seq expression matrix. DropoutKiller does **not** estimate a fully imputed matrix. It estimates a sparse event set

`D = {(g,c): X_gc = 0 and evidence supports technical dropout}`

and only those coordinates are eligible for replacement.

The biological-prior branch is restricted to **PPI and pathway information**. GRN / TF-target priors are intentionally excluded from the package contract.

## 2. SuperCell-style membership

For a low-dimensional cell representation `z_c`, memberships are constructed separately inside hard strata such as major cell type and, optionally, condition or donor.

Within stratum `s`, construct a Euclidean kNN graph `G_s=(V_s,E_s)` and apply walktrap hierarchical community detection. With `n_s` cells and graining level `gamma`, the requested number of memberships is

`K_s = max(1, round(n_s/gamma))`.

The dendrogram is cut at `K_s`, or at least the number of disconnected graph components. For large strata an anchor subset is clustered first; omitted cells are assigned to the nearest membership centroid in embedding space.

This follows the computational structure of SuperCell but moves biological-stratum separation **before** graph construction so that cross-stratum edges cannot influence membership boundaries.

## 3. Membership-local low-rank model

For membership `k`, let

`Y_k = X[, C_k]`.

Compute an uncentered rank-`r_k` approximation

`Y_k ~= U_k Sigma_k V_k^T = Yhat_k`.

When genes greatly outnumber cells, the implementation may use the mathematically equivalent Gram decomposition

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

Observed non-zero expression is therefore excluded before any recovery prior is consulted.

## 5. Confidence against the biological-zero null

Negative low-rank values estimate the symmetric reconstruction-error distribution associated with biological zeros. The working null is

`E_gk ~ N(0, sigma_gk^2)`

with

`sigma_gk^2 = mean(e^2 : e = Yhat_gj < 0)`.

For an ALRA-gated zero,

`C_gc = Phi(Yhat_gc / sigma_gk)`.

If too few negative values are available to estimate `sigma_gk`, the diagnostic candidate is assigned neutral confidence `0.5`; it cannot enter the default `0.95` high-confidence mask.

`C_gc` is a null-tail confidence score, **not** a calibrated posterior `P(dropout | data)`.

The default mask is

`M_gc = I(X_gc=0) I(Yhat_gc>tau_gk) I(C_gc>=0.95)`.

## 6. Cell-space weighted borrowing

For masked event `(g,c)`, only cells in the same membership are eligible donors. With squared latent distance

`d_cj^2 = ||z_c-z_j||_2^2`,

Gaussian weights are

`w_cj proportional to exp(-d_cj^2 / sigma_c^2)`.

Neighbor lookup is performed only for cells carrying masked events. The default adaptive `sigma_c` is the median positive distance among selected neighbors.

By default, donors with `X_gj=0` are excluded and weights are renormalized:

`Xhat_gc^cell = sum_j w_cj X_gj I(X_gj>0) / sum_j w_cj I(X_gj>0)`.

Setting `neighbor_positive_only=FALSE` restores the literal all-neighbor Gaussian average.

## 7. PPI/pathway biological prior

Two optional sources are supported:

- `A^PPI`: PPI adjacency.
- `A^path`: pathway-derived gene interaction adjacency.

No GRN / TF-target prior is accepted.

Each supplied prior is aligned to expression genes, self-edges are removed, and every target row is L1-normalized:

`sum_h |A_gh| = 1` for rows with support.

When both PPI and pathway priors are supplied, they are fused as

`A* = omega_PPI A^PPI + omega_path A^path`

with

`omega_PPI >= 0`, `omega_path >= 0`, and `omega_PPI + omega_path = 1`.

The fused matrix is row-normalized again. If weights are not supplied, available priors receive equal weight.

Direct raw-expression graph averaging is scale-sensitive. Therefore, within membership `k`:

`Z_hc = (X_hc - mu_hk) / s_hk`.

The prior prediction is

`Zhat_gc = sum_h A*_gh Z_hc`

and is mapped back to the target-gene scale:

`Xhat_gc^prior = max(0, mu_gk + s_gk Zhat_gc)`.

The sparse implementation algebraically rewrites the operation so that a global dense standardized `G x C` matrix is not materialized.

The prior branch is **downstream of the dropout gate**. PPI/pathway information cannot convert an otherwise unsupported biological zero into a dropout candidate.

## 8. Selective hybrid recovery

For masked events where both predictions exist,

`Xhat_gc = alpha Xhat_gc^cell + (1-alpha) Xhat_gc^prior`

with default `alpha = 0.75`.

Thus the default recovery is

`75% cell borrowing + 25% PPI/pathway prior`.

If only one component exists, its weight is renormalized to 1. If neither exists, the zero remains unchanged.

Finally,

`Xfinal_gc = X_gc` for `M_gc = 0`

and

`Xfinal_gc = Xhat_gc` for `M_gc = 1` when a positive prediction is available.

## 9. Invariants

1. Mask entries must correspond to original zeros.
2. Observed non-zero values are unchanged exactly.
3. Cell-space recovery never crosses membership boundaries.
4. Supplied hard biological strata cannot be crossed during membership graph construction.
5. PPI/pathway priors are optional and never participate in dropout detection.
6. GRN/TF-target priors are not part of the public or internal biological-prior contract.
7. Recovered values are continuous expression estimates, not raw counts.
8. Sparse input remains sparse in the end-to-end recovery path.
9. Recovery is one-pass; recovered values are never fed back into detection or membership construction.

## 10. Statistical trade-offs

- Hard strata prevent cross-label borrowing but make recovery dependent on annotation quality.
- The Gaussian negative-tail model is a working null; its confidence is not an FDR or posterior probability.
- Positive-only cell borrowing targets the conditional-positive local mean and may be upward-biased.
- PPI and pathway resources may include relationships that are valid globally but inactive in the assayed cell state. This is why the prior branch is optional, applied only after the expression-derived dropout gate, and receives only `1-alpha` weight when cell borrowing is available.
- Pathway edge direction is source-resource dependent. The edge-list constructor defaults pathway edges to directed and PPI edges to symmetric, but users should set `directed` according to the semantics of the resource.
- The pipeline deliberately avoids iterative re-imputation because iterative reuse can amplify the package's own predictions.
