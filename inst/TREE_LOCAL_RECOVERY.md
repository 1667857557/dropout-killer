# Tree-aware local recovery contract

## Statistical target

For a selected technical-dropout event `(g,c)`, version 0.6 targets

\[
E[X_{gc}\mid X_{gc}>0,\;\mathcal T,\;z_c,\;X_{c,-g},\;R,\;m_c,\;s_c],
\]

where `T` is the retained SuperCell walktrap hierarchy, `z_c` is the biological embedding, `R` is the reliable-coordinate mask, `m_c` is the final SuperCell membership, and `s_c` is an optional harder biological stratum.

The final membership is the recovery borrowing block. Tree and embedding distances determine continuous donor importance **inside** that block.

## Information hierarchy

```text
hard biological stratum
        >
final SuperCell membership
        >
walktrap hierarchy within membership
        >
continuous embedding distance
        >
positive target support
        >
non-target coexpression residual
```

A donor outside the final membership has exactly zero recovery weight. An explicit hard biological stratum may further split a membership, but never expands borrowing beyond the membership.

## Tree distance

For two hierarchy leaves `c` and `j`, let `LCA(c,j)` be their lowest common ancestor. The implementation indexes the walktrap merge tree and assigns each internal node a scale derived from its subtree size:

\[
d_T(c,j)=\frac{\log |\mathrm{subtree}(LCA(c,j))|}{\log n_s}\in[0,1].
\]

Small early-merged subtrees therefore have small distance. Approximate SuperCell builds retain the anchor hierarchy; pairs without exact tree leaves use the embedding term only.

## Membership-local donor kernel

For cells in the same final recovery block:

\[
w_{cj}=\exp\left[
-\alpha\frac{d_T(c,j)}{\tau_T}
-(1-\alpha)\frac{d_E(c,j)^2}{2h_c^2}
\right],
\]

where `alpha = tree_weight`, `tau_T = tree_tau`, and `h_c` is the distance to the `local_k`-th candidate neighbor inside the membership.

The candidate set is restricted to the final membership and optionally capped by `candidate_k`. The sparse weight matrix is built once per recovery run. Its elementwise square, row-weight sums, adaptive bandwidths, and cell-level weighted tree/embedding diagnostics are cached once and reused across genes.

## Query-specific positive baseline

For reliable positive donors inside the membership,

\[
D_{gc}=\{j:m_j=m_c,\;X_{gj}>0,\;R_{gj}=1\},
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

Local positive prevalence is

\[
\pi_{gc}^{local}=\frac{\sum_jw_{cj}I(X_{gj}>0)}{\sum_jw_{cj}}.
\]

For engineering efficiency, these quantities are computed for target genes in batches using matrix multiplication:

\[
N=XW^T,\qquad D=BW^T,\qquad S_2=X^{\circ2}W^T,\qquad D_2=BW^{\circ2T},
\]

where `B=I(X>0)`. Then

\[
\mu=N\oslash D,
\qquad
n_{eff}=D^{\circ2}\oslash D_2.
\]

This removes per-event local-statistic loops.

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

## One residual model per gene and membership

The original tree-local implementation fit a separate weighted ridge for every dropout event. With millions of high-confidence events, runtime therefore scaled with event count even though cells were already partitioned into small memberships.

The batched implementation fits one residual model per target gene and recovery block.

For target gene `g`, let

\[
Q_g=\{c:(g,c)\text{ is a selected dropout event in the membership}\}.
\]

Donor importance is aggregated over the actual target query cells:

\[
\bar w_{gj}=\sum_{c\in Q_g}w_{cj}.
\]

This retains the principle that donors close to the dropout queries contribute more, but avoids fitting the same gene separately for every query cell.

The residual ridge is

\[
\widehat\beta_g=(F^T\bar W_gF+P)^{-1}F^T\bar W_gr_g,
\]

where the intercept is unpenalized and `P` contains the ridge penalty on factor coefficients.

The weighted linear smoother provides analytic leave-one-out residual predictions. A gene-level held-out factor coefficient is

\[
q_g^{pred}=\mathrm{clip}_{[0,1]}
\frac{\sum_j\bar w_{gj}\widehat r_{gj}^{(-j)}r_{gj}}
{\sum_j\bar w_{gj}(\widehat r_{gj}^{(-j)})^2}.
\]

Query-specific local information then supplies

\[
q_{gc}^{info}=\frac{n_{eff,gc}}{n_{eff,gc}+\kappa},
\]

so

\[
q_{gc}=q_g^{pred}q_{gc}^{info}.
\]

The final recovery is

\[
\widehat X_{gc}=\max\{0,\mu_{gc}^{local,+}+q_{gc}\widehat r_{gc}\}.
\]

If coexpression has no held-out gain, the model collapses to the query-specific local positive mean rather than the whole-membership mean.

## Predictive uncertainty

For a given gene-level residual model, the deployed shrinkage remains query-specific because `n_eff` varies by query. Therefore held-out residual SSE is evaluated for each query's actual `q_gc` using its local donor weights without refitting the ridge.

Local-mean estimation uncertainty scales with `1/n_eff`, and factor uncertainty uses the gene-level weighted ridge leverage. The working predictive variance is

\[
V_{gc}\approx\widehat\sigma_{g,LOO,c}^2(1+q_{gc}^2h_{gc})+\frac{s_{gc,local}^2}{n_{eff,gc}}.
\]

This remains an approximate predictive variance and must be checked by held-out coverage.

Positive multiple-imputation draws keep the existing Gamma moment match:

\[
X_{gc}^{(b)}\sim Gamma\left(\frac{m_{gc}^2}{V_{gc}},\frac{V_{gc}}{m_{gc}}\right),
\]

where the second parameter is scale. Therefore the sampled mean and variance equal the stored predictive moments while remaining positive.

## Computational contract

For each final membership block:

```text
build W once
cache W^2 and cell-level geometry summaries
compute factor scores once
batch target-gene local statistics
for each target gene:
    fit one residual ridge
    predict every dropout query for that gene
index predictions back into the event table
```

If `E` is the number of selected dropout coordinates, `G_t` is the number of target genes, `n_m` is membership size, and `K` is factor rank, the expensive regression work changes from approximately

\[
O(E\,n_mK^2)
\]

to approximately

\[
O(G_t\,n_mK^2),
\]

plus batched matrix multiplications and final event indexing. Event count no longer determines the number of ridge solves.

## Required benchmark before default promotion

`tree_local_factor` remains a benchmark engine until independent validation establishes superiority. The same oracle masks must compare:

1. positive membership mean;
2. embedding-only kernel mean;
3. tree-only weighted mean;
4. tree + embedding weighted mean;
5. current masked factor;
6. tree-local mean + gene-level residual factor.

Benchmarks must include random observed-positive masking, original-UMI count strata, Binomial-zero, and full binomial thinning followed by re-normalization/PCA/SuperCell reconstruction. Report RMSE, MAE, bias, Pearson/CCC, prediction-interval coverage, gene-variance error, distributional distance, tree-distance strata, effective-donor strata, and biological strata.

The engine should become default only if it improves held-out recovery while preserving bias, calibration, differential variability, and hard biological boundaries across repeated random seeds.