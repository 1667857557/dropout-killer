# Tree-aware local recovery contract

## Statistical target

For a selected technical-dropout event `(g,c)`, version 0.6 adds a recovery engine targeting

\[
E[X_{gc}\mid X_{gc}>0,\;\mathcal T,\;z_c,\;X_{c,-g},\;R,\;s_c],
\]

where `T` is the retained SuperCell walktrap hierarchy, `z_c` is the biological embedding, `R` is the reliable-coordinate mask, and `s_c` is a hard biological stratum. The final gamma membership is treated as a high-weight tree cut, not as an assumption that every cell in that membership is exchangeable.

## Information hierarchy

```text
hard biological stratum
        >
walktrap hierarchy
        >
final gamma cut
        >
continuous embedding distance
        >
positive target support
        >
non-target coexpression residual
```

Different hard strata have exactly zero donor weight. Within a stratum, nearby sibling memberships can contribute when they are close in the retained hierarchy/embedding.

## Tree distance

For two hierarchy leaves `c` and `j`, let `LCA(c,j)` be their lowest common ancestor. The implementation indexes the walktrap merge tree and assigns each internal node a scale derived from its subtree size:

\[
d_T(c,j)=\frac{\log |\mathrm{subtree}(LCA(c,j))|}{\log n_s}\in[0,1].
\]

Small early-merged subtrees therefore have small distance. Cells in disconnected hierarchy components receive maximal tree distance. Approximate SuperCell builds retain the anchor hierarchy; pairs without exact tree leaves use the embedding term only.

## Joint donor kernel

For cells in the same hard stratum:

\[
w_{cj}=\exp\left[
-\alpha\frac{d_T(c,j)}{\tau_T}
-(1-\alpha)\frac{d_E(c,j)^2}{2h_c^2}
\right],
\]

where `alpha = tree_weight`, `tau_T = tree_tau`, and `h_c` is the distance to the `local_k`-th candidate neighbor. If tree distance is unavailable for a pair, its effective tree weight is zero instead of fabricating a hierarchy distance.

Candidate donors are the union of the query's final membership and embedding-nearest neighbors inside the hard stratum. The sparse weight matrix is built once and reused across genes.

## Query-specific positive baseline

For reliable positive donors

\[
D_{gc}=\{j:X_{gj}>0,R_{gj}=1,s_j=s_c\},
\]

the local baseline is

\[
\mu_{gc}^{local,+}=\frac{\sum_{j\in D_{gc}}w_{cj}X_{gj}}{\sum_{j\in D_{gc}}w_{cj}}.
\]

Unlike a positive membership mean, this baseline varies across query cells inside the same membership.

The effective donor size is

\[
n_{eff,gc}=\frac{(\sum_jw_{cj})^2}{\sum_jw_{cj}^2}.
\]

The weighted local variance is

\[
s_{gc}^2=\frac{\sum_jw_{cj}(X_{gj}-\mu_{gc})^2}{W_1-W_2/W_1},
\]

with `W1=sum(w)` and `W2=sum(w^2)`.

Local positive prevalence is also stored:

\[
\pi_{gc}^{local}=\frac{\sum_jw_{cj}I(X_{gj}>0)}{\sum_jw_{cj}}.
\]

## Leave-one-out local residuals

The sparse kernel has zero diagonal, so a reliable positive donor `j` receives a baseline computed without itself:

\[
r_{gj}=X_{gj}-\mu_{gj}^{(-j)}.
\]

All current target genes are excluded from factor-state learning. Thus the prediction path remains

\[
X_{c,-T}\rightarrow f_c\rightarrow\widehat r_{gc},
\]

never `X_gc -> factor -> Xhat_gc`.

## Query-weighted residual ridge

For a query cell, positive residual donors are fit with

\[
\widehat\beta_{gc}=(F^TW_cF+P)^{-1}F^TW_cr_g,
\]

where the intercept is unpenalized. Because `W_c` depends on the query cell, biologically closer donors have more influence on the residual model as well as on the positive baseline.

The weighted linear smoother permits analytic leave-one-out residual predictions. Relative to a zero-residual null, factor contribution is shrunk by the weighted held-out squared-error optimum

\[
q_{gc}^{pred}=\mathrm{clip}_{[0,1]}
\frac{\sum_jw_{cj}\widehat r_{gj}^{(-j)}r_{gj}}
{\sum_jw_{cj}(\widehat r_{gj}^{(-j)})^2}.
\]

Local information contributes an additional shrinkage

\[
q_{gc}^{info}=\frac{n_{eff,gc}}{n_{eff,gc}+\kappa},
\]

so

\[
q_{gc}=q_{gc}^{pred}q_{gc}^{info}.
\]

The final recovery is

\[
\widehat X_{gc}=\max\{0,\mu_{gc}^{local,+}+q_{gc}\widehat r_{gc}\}.
\]

If coexpression has no local held-out gain, the model collapses to the query-specific local positive mean rather than the whole-membership mean.

## Predictive uncertainty

The local residual variance is estimated from weighted held-out residual error. Local-mean estimation uncertainty scales with `1/n_eff`, and query factor uncertainty uses weighted ridge leverage. The working predictive variance is

\[
V_{gc}\approx\widehat\sigma_{g,LOO}^2(1+q_{gc}^2h_{gc})+\frac{s_{gc,local}^2}{n_{eff,gc}}.
\]

This remains an approximate predictive variance and must be checked by held-out coverage.

Positive multiple-imputation draws keep the existing Gamma moment match:

\[
X_{gc}^{(b)}\sim Gamma\left(\frac{m_{gc}^2}{V_{gc}},\frac{V_{gc}}{m_{gc}}\right),
\]

where the second parameter is scale. Therefore the sampled mean and variance equal the stored predictive moments while remaining positive.

## Required benchmark before default promotion

`tree_local_factor` is introduced as a parallel engine, not immediately promoted over `masked_factor`. The same oracle masks must compare:

1. positive membership mean;
2. embedding-only kernel mean;
3. tree-only weighted mean;
4. tree + embedding weighted mean;
5. current masked factor;
6. local mean + residual factor;
7. tree-local weighted residual factor.

Benchmarks must include MCAR positive masking, original-UMI count strata, Binomial-zero, and full binomial thinning followed by re-normalization/PCA/SuperCell reconstruction. Report RMSE, MAE, bias, Pearson/CCC, prediction-interval coverage, gene-variance error, distributional distance, tree-distance strata, effective-donor strata, and B-cell biological strata. A post-recovery rise in coexpression is not independent validation because coexpression participates in prediction.

The engine should become default only if it improves held-out recovery while preserving bias, calibration, differential variability, and hard biological boundaries across repeated random seeds.
