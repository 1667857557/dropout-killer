# DropoutKiller 0.8.1 mathematical contract

## 1. Estimand

For raw UMI counts \(X\in\mathbb N_0^{G\times n}\), DropoutKiller separates two questions:

1. **Detection:** which observed RNA zeros have sufficient cell-specific evidence for latent expression?
2. **Recovery:** conditional on a zero being selected, what positive expression value is supported by target-safe cell-state information?

The method does not assume that every observed zero is a separate technical-dropout process, and it does not use scGACL or a Gamma-Normal gene-by-cell-type dropout posterior.

Observed nonzero coordinates are immutable.

---

## 2. Working scale

The production workflow transforms raw counts once:

\[
Y_{gc}=\log\left(1+s\frac{X_{gc}}{\sum_hX_{hc}}\right),\qquad s=10^4.
\]

All detection and deterministic recovery values live on this normalized log-expression scale.

---

## 3. Broad biological boundaries

The user supplies a broad biological class \(B_c\) for each cell. These classes are hard borrowing boundaries. Fine subclusters are not inferred as boundaries.

Within each class, a kNN graph is constructed and Walktrap retains its complete merge history. The target number of final memberships is

\[
K_B=\max\left(1,\left\lfloor\frac{n_B}{\gamma}\right\rfloor,C_B\right),
\]

where \(C_B\) is the connected-component count and the default is \(\gamma=150\).

No graph edge or recovery donor is introduced across different broad classes.

---

## 4. RNA low-rank evidence

For normalized expression, a truncated SVD provides

\[
Y\approx U_kD_kV_k^T=\widehat Y.
\]

Automatic rank uses the retained ALRA singular-value-spacing heuristic. For each gene, the reconstructed negative tail provides a scale estimate. The detector feature \(z_{gc}\) is the standardized margin of \(\widehat Y_{gc}\) above the gene-specific ALRA-style gate.

This quantity is latent-expression evidence. It is not a posterior probability.

---

## 5. Hierarchy-local support

For query cell \(c\), donors are selected from its own final membership first. Only when fewer than `neighbor_k` usable donors exist are cells borrowed from the nearest connected memberships in the retained Walktrap hierarchy.

Let \(W_{cj}\) be the row-normalized donor weight and \(I(Y_{gj}>0)\) indicate detected positive expression. Then

\[
s_{gc}=\sum_jW_{cj}I(Y_{gj}>0)
\]

and the local support feature is

\[
L_{gc}=\operatorname{logit}\left(\frac{s_{gc}+0.01}{1.02}\right).
\]

Self-neighbors are excluded.

---

## 6. Membership context

For query membership \(m\) inside broad group \(B\), define \(k_m\) and \(n_m\) as positive and total cell counts for the gene within the membership and \(k_B,n_B\) analogously for the broad group. The context feature is

\[
H_{gc}=\frac{n_m}{n_m+50}
\left[
\operatorname{logit}\left(\frac{k_m+0.5}{n_m+1}\right)-
\operatorname{logit}\left(\frac{k_B+0.5}{n_B+1}\right)
\right].
\]

The multiplicative term shrinks unstable small-membership contrasts toward zero.

---

## 7. Calibrated Lean detector

The production detector is a class-balanced logistic model on exactly three features:

\[
\eta_{gc}=\beta_0+\beta_1z_{gc}+\beta_2L_{gc}+\beta_3H_{gc}.
\]

A zero is selected when

\[
\eta_{gc}>t,
\]

where \(t\) is fitted from artificial-thinning calibration. The stored sigmoid \(\sigma(\eta_{gc})\) is a score only; it is not interpreted as a biological dropout posterior.

### Calibration labels

Counts are independently thinned:

\[
X'_{gc}\mid X_{gc}\sim\operatorname{Binomial}(X_{gc},q),
\]

with default \(q=0.5\) and seeds 10001:10003.

- Artificially lost positives: \(X_{gc}>0\) and \(X'_{gc}=0\).
- Proxy negatives: original zeros with prevalence \(\le0.005\) in the current broad class and \(\ge0.20\) in another broad class.

Every thinning replicate rebuilds the detector geometry. The default `fpr=0.01` controls the empirical training proxy-negative exceedance rate only. It is neither a biological FDR nor a guarantee on independent data.

---

## 8. Paired RNA+ATAC geometry

For paired multiome Seurat inputs, RNA PCA and ATAC TF-IDF/LSI are rebuilt. LSI component 1 is excluded. WNN is recomputed separately inside each broad biological class.

The affinity follows the SuperCell 2.0 multimodal-neighbor convention

\[
A_{ij}=1-2d_{ij}^2,
\]

followed by symmetrization and removal of self-loops. Walktrap uses this affinity directly. Missing WNN edges are not replaced with RNA-only Euclidean edges.

ATAC therefore informs neighborhood geometry, not the recovered RNA magnitude directly.

---

## 9. Selective P1 stabilized-state recovery

Let \(M\) be the detector mask. Recovery is evaluated only where \(M_{gc}=1\).

Target genes are deterministically partitioned into folds. For a target fold \(T_f\), every gene in \(T_f\) is excluded from predictor-state construction, preventing the path

\[
Y_g\rightarrow\text{cell state}\rightarrow\widehat Y_g.
\]

For non-target genes, standardized predictor expression receives one row-stochastic geometry smoothing step

\[
Z^*=(1-\rho)Z+\rho ZP^T,
\]

with default \(\rho=0.25\), after which a low-dimensional factor state is learned.

For target gene \(g\), only positive donor cells are used. With factor design matrix \(Q_g\), the ridge state is

\[
\widehat\beta_g=(Q_g^TQ_g+P_\lambda)^{-1}Q_g^Ty_g,
\]

with an unpenalized intercept and default ridge penalty \(\lambda=2\).

Analytic leave-one-out predictions estimate whether the factor state improves over the positive-donor mean. Unsupported factor contributions shrink toward the conservative mean state. The production defaults use five target-gene folds, support-adaptive factor rank, and `bias_kappa=10`.

---

## 10. Output invariants

If \(Y^{out}\) is the deterministic recovered matrix,

\[
Y^{out}_{gc}=Y_{gc}\qquad\text{for every }Y_{gc}>0.
\]

Only selected observed zeros may change. Recovery values are constrained non-negative on the working scale. Predictive variance is stored for recovery engines that provide an uncertainty model.

---

## 11. Validation contract

Detector fitting and performance evaluation must use different masks or datasets. A valid benchmark should report separately:

- recoverable-zero precision/recall or analogous detection metrics;
- reconstruction MSE/MAE/Spearman conditional on held-out positives;
- stability of DE and gene-set enrichment;
- clustering and marker preservation;
- uncertainty coverage when predictive variance is used.

An increase in post-recovery coexpression is not self-validating because coexpression contributes to the recovery state.

---

## 12. Supported production API

The production detector is only the SuperCell hierarchy Lean detector, with RNA and WNN modality variants. Historical scGACL, empirical-Bayes zero-null, local ALRA-quantile, and global-ALRA detector routes are not part of the current API. Comparator recovery engines may remain internally available for controlled ablation, but the production recovery default is P1 stabilized state.
