# DropoutKiller mathematical and engineering contract

## 1. Scope

Let `X in R_+^(G x C)` be a normalized scRNA-seq expression matrix. DropoutKiller does **not** estimate a fully imputed matrix. It estimates a sparse event set

`D = {(g,c): X_gc = 0 and evidence supports technical dropout}`

and only those coordinates are eligible for replacement.

## 2. SuperCell-style membership

For a low-dimensional cell representation `z_c`, memberships are constructed separately inside hard strata such as major cell type (and optionally condition/donor).

Within stratum `s`, construct a Euclidean kNN graph `G_s=(V_s,E_s)` and apply walktrap hierarchical community detection. With `n_s` cells and graining level `gamma`, the requested number of memberships is

`K_s = max(1, round(n_s/gamma))`.

The dendrogram is cut at `K_s` (or at least the number of disconnected graph components). For large strata an anchor subset is clustered first; omitted cells are assigned to the nearest membership centroid in embedding space.

This follows the computational structure of SuperCell but deliberately moves biological-stratum separation **before** graph construction instead of relying only on post-hoc splitting.

## 3. Membership-local low-rank model

For membership `k`, let

`Y_k = X[, C_k]`.

Compute an uncentered rank-`r_k` SVD approximation

`Y_k ~= U_k Sigma_k V_k^T = Yhat_k`.

A numeric rank can be supplied. With `rank="auto"`, singular-value spacings are compared with their tail distribution, following the ALRA rank-selection principle; the last spacing more than `rank_z` tail standard deviations above noise is retained.

## 4. ALRA-derived zero-preserving gate

For each gene `g` in membership `k`, define

`q_gk = Q_p(Yhat_g,Ck)`

and, only when the lower quantile is negative,

`tau_gk = |q_gk|`.

An observed zero becomes an ALRA-gated candidate only if

`X_gc = 0` and `Yhat_gc > tau_gk`.

This is stricter than the original prototype, which applied a sigmoid to every matrix entry and therefore could mark observed non-zero values.

## 5. Confidence against the biological-zero null

ALRA motivates using negative low-rank values to estimate the symmetric reconstruction-error distribution associated with biological zeros. DropoutKiller uses the working null

`E_gk ~ N(0, sigma_gk^2)`

with

`sigma_gk^2 = mean(e^2 : e = Yhat_gj < 0)`.

For an ALRA-gated zero,

`C_gc = Phi(Yhat_gc / sigma_gk)`.

If there are too few negative values to estimate `sigma_gk`, the event is retained in the diagnostic candidate table but assigned neutral confidence `0.5`, so it does not enter the default high-confidence mask. `C_gc` is a null-tail confidence score, **not a calibrated posterior P(dropout | data)**. The default high-confidence mask is

`M_gc = I(X_gc=0) I(Yhat_gc>tau_gk) I(C_gc>=0.95)`.

## 6. Cell-space weighted borrowing

For candidate `(g,c)`, consider only cells in the same membership. With squared latent distance

`d_cj^2 = ||z_c-z_j||_2^2`, 

Gaussian weights are

`w_cj proportional to exp(-d_cj^2 / sigma_c^2)`.

Neighbor lookup is performed only for cells carrying masked events (RANN exact kNN within the membership), avoiding an O(n_k^2) full distance matrix. The default adaptive `sigma_c` is the median positive distance among the selected neighbors. For recovery of gene `g`, donors with `X_gj=0` are excluded by default (`neighbor_positive_only=TRUE`). This is configurable; setting it to `FALSE` recovers the literal all-neighbor Gaussian average. The default renormalized prediction is

`Xhat_gc^cell = sum_j w_cj X_gj I(X_gj>0) / sum_j w_cj I(X_gj>0)`.

This avoids reinforcing a technical zero by averaging it with other unresolved zeros.

## 7. Gene-network prior

Let `A` be a PPI/pathway/GRN adjacency matrix aligned to expression genes. Self-edges are removed and each target row is L1-normalized.

Direct raw-expression graph averaging is scale-sensitive. Therefore, inside each membership, genes are standardized:

`Z_hc = (X_hc-mu_hk)/s_hk`.

The network prediction is

`Zhat_gc = sum_h A_gh Z_hc`,

then mapped back to the target-gene scale:

`Xhat_gc^gene = max(0, mu_gk + s_gk Zhat_gc)`.

The implementation algebraically rewrites this operation so sparse adjacency and sparse expression can be multiplied without materializing a full standardized `G x C` matrix.

## 8. Selective hybrid recovery

For mask entries, when both predictions exist,

`Xhat_gc = alpha Xhat_gc^cell + (1-alpha) Xhat_gc^gene`,

with default `alpha=0.75`.

If only one component is available, component weights are renormalized rather than silently replacing the unavailable component with zero. If neither component is available, the observed zero remains zero.

Finally,

`Xfinal_gc = X_gc` for `M_gc=0`,

and

`Xfinal_gc = Xhat_gc` for `M_gc=1` when a positive prediction is available.

## 9. Invariants

The package enforces the following invariants:

1. Mask entries must correspond to original zeros.
2. Observed non-zero values are unchanged exactly by the core workflow.
3. Recovery is membership-local in cell space.
4. Supplied hard strata cannot be crossed during membership graph construction.
5. Gene priors are optional and never required for cell-space recovery.
6. Recovered values are continuous expression estimates and are not relabeled as raw counts.
7. Sparse input remains sparse in the end-to-end recovery path; dense `G x C` prediction matrices are not created globally.


## 10. Statistical trade-offs

- Hard biological strata prevent cross-label contamination but make membership quality dependent on annotation quality.
- The Gaussian negative-tail model is a working null. Its confidence score is not an FDR or posterior probability.
- Positive-only cell borrowing protects against unresolved technical zeros but targets the conditional-positive local mean and may be upward-biased.
- Gene networks may encode context-mismatched relationships. For that reason the network branch is optional, follows the dropout gate rather than determining it, and receives only `1-alpha` weight when both components are available.
- The pipeline is deliberately one-pass: recovered values never enter a second detection or membership round.
