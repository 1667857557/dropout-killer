# DropoutKiller Benchmark Design

## 1. Purpose

This document defines the benchmark used to decide whether a recovery engine actually reconstructs expression that is closer to the original biological signal.

The central benchmark is deliberately simple and falsifiable:

> Start from expression values that were genuinely observed, randomly hide a subset of positive gene-by-cell coordinates by setting them to zero, give the recovery engine the exact artificial-dropout mask, and compare the recovered values with the untouched original values.

This creates an **oracle recovery benchmark** because the ground-truth expression at every artificial dropout is known exactly.

The benchmark must separate two questions:

1. **Recovery:** if the location of a technical dropout is already known, how accurately can its expression magnitude be reconstructed?
2. **Detection + recovery:** after sequencing loss creates new zeros, can the complete pipeline detect the correct zeros and recover their expression?

These two questions must not be conflated. The primary benchmark therefore freezes the original biological geometry and supplies the true artificial mask. Count thinning is a separate secondary benchmark.

---

# 2. Primary benchmark: random observed-positive knockout

## 2.1 Reference matrix

Let

\[
X\in\mathbb R_+^{G\times C}
\]

be the normalized/log-transformed expression matrix used by DropoutKiller, with genes in rows and cells in columns.

The untouched matrix is the reference:

\[
X^{ref}=X.
\]

Only coordinates with observed positive expression are eligible for artificial masking:

\[
\mathcal E
=
\{(g,c):X^{ref}_{gc}>0\}.
\]

This restriction is essential. A naturally observed zero has no known non-zero ground truth, so it cannot be used to score recovery magnitude.

---

## 2.2 Artificial dropout mask

For an eligible coordinate \((g,c)\), draw

\[
M_{gc}\sim Bernoulli(p_{mask}).
\]

Then construct

\[
X^{mask}_{gc}
=
\begin{cases}
0,&M_{gc}=1\\
X^{ref}_{gc},&M_{gc}=0.
\end{cases}
\]

The exact artificial-dropout truth is retained separately:

\[
T_{gc}=X^{ref}_{gc},\qquad M_{gc}=1.
\]

The recovery algorithm receives:

- `x = X_mask`;
- the exact sparse logical `mask = M`;
- the original biological embedding;
- the original SuperCell membership/hierarchy;
- the same hard biological strata used in the real analysis.

It must **not** receive the hidden values \(T_{gc}\).

---

# 3. Why the embedding and hierarchy are frozen in the primary benchmark

The primary benchmark is intended to isolate recovery error.

Therefore construct the biological state representation from the untouched data once:

```text
X_ref
  |
  +--> PCA / Harmony / WNN representation Z_ref
  |
  +--> SuperCell kNN graph
  |
  +--> walktrap hierarchy T_ref
  |
  +--> final membership cut M_ref
```

After artificial masking:

```text
X_mask + Z_ref + T_ref + M_ref + oracle dropout mask
                         |
                         v
                      recovery
                         |
                         v
                       X_hat
```

Do **not** recompute PCA, membership, or the hierarchy after the random knockout in this primary experiment.

Otherwise the measured error becomes a mixture of:

\[
\text{recovery error}
+
\text{embedding perturbation}
+
\text{graph perturbation}
+
\text{membership perturbation}.
\]

That is useful for an end-to-end stress test, but not for determining whether the recovery estimator itself is mathematically better.

---

# 4. Masking design must reflect single-cell biology

Uniformly sampling all positive coordinates is useful but insufficient because scRNA-seq expression is strongly heterogeneous across genes, cells, cell states, and abundance levels.

The benchmark must therefore report both an overall result and stratified results.

## 4.1 Recommended mask fractions

Use several independent difficulty levels:

\[
p_{mask}\in\{0.05,0.10,0.20\}.
\]

The 5% mask is a low-perturbation benchmark. The 10% mask is the primary default. The 20% mask is a stress test.

A benchmark run must never silently mask an excessive fraction of one sparse gene or one cell.

Recommended constraints:

```text
maximum masked fraction per gene <= 20% of its observed-positive cells
maximum masked fraction per cell <= 20% of its observed-positive genes
```

If a requested random draw violates a constraint, resample the offending coordinates.

---

## 4.2 Positive-support strata

A recovery method has fundamentally different information when a target gene is positive in 5 nearby cells versus 100 nearby cells.

For each artificial event calculate the number or Kish effective number of eligible positive donors and stratify the results, for example:

```text
n_eff < 3
3 <= n_eff < 5
5 <= n_eff < 10
10 <= n_eff < 20
n_eff >= 20
```

Unsupported events should not be discarded from the benchmark. They must be reported as their own stratum because failure to recover sparse genes is part of the method's operating characteristics.

---

## 4.3 Expression-magnitude strata

Observed positives should be divided by their original reference magnitude, preferably using quantiles within gene or within biological stratum.

A simple reporting scheme is:

```text
Q1: low positive expression
Q2
Q3
Q4: high positive expression
```

If raw UMI counts are available, additionally report:

```text
count = 1
count = 2
count = 3-5
count > 5
```

This is biologically important because sequencing-induced zeros are disproportionately generated from low-abundance transcripts.

---

## 4.4 Gene detection-rate strata

For gene \(g\), define the positive detection rate inside the relevant hard biological stratum:

\[
\pi_g
=
\frac{1}{C_s}
\sum_{c\in s}I(X_{gc}>0).
\]

Suggested bins:

```text
5-10%
10-20%
20-50%
50-80%
>80%
```

Genes below the minimum support required by a recovery engine may be retained in a separate `insufficient_support` category rather than being silently removed.

---

## 4.5 Biological-state strata

Report results separately for the available trusted cell states, for example in a B-cell dataset:

```text
Naive B
Memory B
Intermediate B
Plasma
```

and also by final SuperCell membership.

This prevents an apparently good global RMSE from hiding over-smoothing of a minority biological state.

---

# 5. Balanced artificial-dropout sampling

A purely coordinate-uniform mask is dominated by genes that are common and highly detected.

The preferred benchmark therefore contains two sampling schemes.

## 5.1 Coordinate-MCAR benchmark

Sample uniformly from all eligible positive coordinates.

This answers:

> What is the average recovery error over the observed-positive matrix?

## 5.2 Stratified balanced benchmark

Sample approximately equal numbers of artificial events across combinations of:

- biological state;
- gene detection-rate bin;
- original expression-magnitude bin;
- positive-donor support bin.

This answers:

> Does the method improve only easy/high-expression genes, or does the improvement generalize across biologically relevant regimes?

Both results should be reported.

---

# 6. Repeated random masks

A single random mask is not sufficient.

Use a fixed list of independent seeds, for example:

```r
benchmark_seeds <- 1:20
```

or preferably 30 seeds for a final report.

All recovery methods must receive **exactly the same mask for a given seed**.

This paired design greatly reduces benchmark noise.

For method \(A\) and method \(B\), compare per-seed differences such as:

\[
\Delta RMSE_s
=
RMSE_{A,s}-RMSE_{B,s}.
\]

Do not compare methods generated from different random masks.

---

# 7. Avoid hyperparameter leakage

Parameters such as

\[
\alpha,\tau_T,k_{local},k_{candidate},\kappa,\lambda
\]

must not be selected on the same artificial-dropout events used for the final performance estimate.

Use two independent mask sets:

```text
Development masks
    -> choose hyperparameters

Locked confirmatory masks
    -> final benchmark only
```

For example:

```text
seeds 1-10   : development/tuning
seeds 101-120: locked confirmatory benchmark
```

The confirmatory seeds should not be inspected while tuning the model.

---

# 8. Recovery methods to compare

Every method must use the same artificial mask and the same reference biological geometry when applicable.

Recommended component ablation:

1. `zero` — leave the artificial dropout as zero;
2. global gene positive mean;
3. final-membership positive mean;
4. `neighbor`;
5. current `masked_factor`;
6. embedding-only local positive mean (`tree_weight = 0`);
7. tree-only local positive mean (`tree_weight = 1`);
8. tree + embedding local positive mean with factor refinement disabled;
9. `tree_local_factor`.

The important comparisons are not only against zero, but against progressively stronger baselines.

In particular:

\[
\text{tree_local_factor}
\quad\text{vs}\quad
\text{tree+embedding local mean}
\]

isolates the incremental value of coexpression residual learning.

Similarly:

\[
\text{tree+embedding local mean}
\quad\text{vs}\quad
\text{embedding-only local mean}
\]

isolates the incremental value of the retained SuperCell hierarchy.

---

# 9. Primary event-level metrics

Let \(\mathcal H\) denote the artificial-dropout coordinates.

## 9.1 RMSE

\[
RMSE
=
\sqrt{
\frac{1}{|\mathcal H|}
\sum_{(g,c)\in\mathcal H}
(\widehat X_{gc}-X^{ref}_{gc})^2
}.
\]

This is the primary magnitude-loss metric.

## 9.2 MAE

\[
MAE
=
\frac{1}{|\mathcal H|}
\sum_{(g,c)\in\mathcal H}
|\widehat X_{gc}-X^{ref}_{gc}|.
\]

MAE is less dominated by a small number of extreme errors.

## 9.3 Bias

\[
Bias
=
\frac{1}{|\mathcal H|}
\sum_{(g,c)\in\mathcal H}
(\widehat X_{gc}-X^{ref}_{gc}).
\]

Positive bias is especially important for a positive-conditional recovery method because it can reveal systematic over-imputation.

## 9.4 Pearson correlation

Measures whether recovered values preserve cell/gene-specific ranking on the held-out coordinates.

## 9.5 Spearman correlation

Provides a ranking metric less sensitive to scale and outliers.

## 9.6 Concordance correlation coefficient

Pearson correlation alone can be high even when predictions have the wrong scale or mean.

Concordance correlation should therefore be included to assess simultaneous agreement in correlation, location, and scale.

---

# 10. Do not use MAPE as the main metric

For low positive expression values,

\[
\frac{|\widehat X-X|}{X}
\]

can explode as \(X\rightarrow0\).

MAPE therefore disproportionately weights the weakest positive values and should not be a primary score.

If a relative-error metric is desired, use it only as a stratified secondary diagnostic with a clearly defined denominator floor.

---

# 11. Gene-level distribution preservation

Recovery should reproduce more than the average event value.

For every sufficiently represented gene, compare the held-out reference distribution with the recovered distribution.

## 11.1 Mean error

\[
\Delta\mu_g
=
\widehat\mu_g-\mu_g^{ref}.
\]

## 11.2 Variance error

\[
E_{var,g}
=
\frac{
|\widehat V_g-V_g^{ref}|
}{V_g^{ref}+\epsilon}.
\]

## 11.3 Wasserstein distance

Use the one-dimensional Wasserstein-1 distance:

\[
W_1(P_g^{rec},P_g^{ref}).
\]

This detects changes in the whole expression distribution rather than only its first two moments.

---

# 12. Differential-variability benchmark

The deterministic recovered matrix contains conditional means and will generally contract variance.

Therefore evaluate both:

1. deterministic point recovery;
2. uncertainty-aware completed draws from `sample_dropout_expression()`.

For each gene compare:

\[
V_g^{ref},
\qquad
V_g^{mask},
\qquad
V_g^{point},
\qquad
V_g^{draw}.
\]

A useful normalized error is:

\[
DVerror_g
=
\frac{|V_g^{method}-V_g^{ref}|}{V_g^{ref}+\epsilon}.
\]

Also evaluate whether the rank ordering of gene variances is preserved.

For methods with predictive uncertainty, use repeated completed matrices, for example:

```r
draws <- sample_dropout_expression(fit, n = 20, seed = 1)
```

and increase to 50 or more draws for the final DV report if computationally feasible.

---

# 13. Predictive uncertainty calibration

For every artificial dropout the true value is known, so the event-level predictive uncertainty can be tested directly.

For a stored predictive mean \(m\) and variance \(v\), the positive-target engine uses the Gamma moment match:

\[
shape=\frac{m^2}{v},
\qquad
scale=\frac{v}{m}.
\]

Use this distribution to form nominal prediction intervals.

Report empirical coverage for at least:

```text
50% interval
80% interval
95% interval
```

For nominal level \(q\),

\[
Coverage(q)
=
\frac{1}{|\mathcal H|}
\sum_{(g,c)\in\mathcal H}
I\{X^{ref}_{gc}\in PI_{gc}(q)\}.
\]

Also report average interval width.

A model should not obtain good coverage merely by returning excessively broad intervals.

Coverage should be stratified by:

- effective donor number;
- original expression magnitude;
- tree distance;
- biological state.

---

# 14. Tree-distance benchmark

For each artificial event calculate the distance to the closest eligible positive donor and/or weighted donor-average tree distance.

Divide events into tree-distance quartiles.

The expected qualitative pattern is:

\[
d_T\uparrow
\Rightarrow
RMSE\uparrow,
\]

and

\[
d_T\uparrow
\Rightarrow
prediction\ uncertainty\uparrow.
\]

If tree distance has no relationship with recoverability, or if tree weighting degrades every distance stratum relative to embedding-only recovery, the retained hierarchy is not providing useful additional information.

---

# 15. Effective-donor benchmark

For query \((g,c)\), define the Kish effective positive-donor size

\[
n_{eff,gc}
=
\frac{(\sum_jw_{cj})^2}{\sum_jw_{cj}^2}.
\]

Report accuracy and uncertainty calibration across `n_eff` bins.

A coherent estimator should usually show:

\[
n_{eff}\uparrow
\Rightarrow
RMSE\downarrow
\]

and

\[
n_{eff}\uparrow
\Rightarrow
prediction\_sd\downarrow,
\]

conditional on comparable expression regimes.

Failure of these relationships is a diagnostic that weights or uncertainty are not behaving as intended.

---

# 16. Biological-boundary stress test

Tree-local recovery is allowed to borrow across final memberships only inside a trusted hard biological stratum.

The benchmark should deliberately sample artificial dropouts from cells close to hard-state boundaries and test whether the method remains unbiased within each state.

For example, report separately:

```text
Naive B near Memory B boundary
Memory B near Naive B boundary
Plasma cells
intermediate states
```

The objective is not to maximize cross-state smoothness. The objective is to reconstruct the original value without erasing state-specific expression.

---

# 17. Structured masking stress tests

Independent random masking is the primary benchmark because its truth is clean and interpretation is simple. Additional structured masks test robustness.

## 17.1 Gene-heavy masking

Randomly select a subset of genes and mask a larger fraction of their positive coordinates.

This tests whether recovery collapses when a target gene has fewer donors.

## 17.2 Cell-heavy masking

Randomly select a subset of cells and mask more positive coordinates in those cells.

This mimics poor-capture cells and tests whether the non-target cell-state representation remains useful.

## 17.3 Local-neighborhood masking

Choose local regions of the embedding/tree and mask target-gene positives in a contiguous neighborhood.

This tests whether the model can recover when the nearest donors are themselves missing and it must broaden the borrowing radius.

These are stress tests and should not replace the coordinate-MCAR primary benchmark.

---

# 18. Secondary benchmark: count-level thinning

Randomly replacing a normalized positive value by zero is excellent for evaluating recovery magnitude because the exact hidden value is known, but it is not a complete sequencing model.

If raw UMI counts \(Y\) are available, add a more biologically realistic thinning experiment.

For retention fraction \(\rho\), generate

\[
Y'_{gc}\sim Binomial(Y_{gc},\rho).
\]

Suggested values:

\[
\rho\in\{0.75,0.50,0.25\}.
\]

A known artificial sequencing-zero event is

\[
Y_{gc}>0
\quad\text{and}\quad
Y'_{gc}=0.
\]

The original count-positive coordinate supplies the event truth.

This experiment naturally creates more zeros among low-count transcripts and therefore better mimics UMI loss.

---

# 19. Full end-to-end thinning benchmark

For the secondary benchmark, rebuild the complete analysis from the thinned counts:

```text
raw counts Y
    |
    v
binomial thinning Y'
    |
    v
normalization
    |
    v
PCA / Harmony / WNN
    |
    v
SuperCell graph + hierarchy + membership
    |
    v
dropout detection
    |
    v
recovery
```

This benchmark measures:

\[
\text{detection}
+
\text{state reconstruction}
+
\text{recovery}
\]

and must therefore report detection and recovery separately.

---

# 20. Detection metrics in the end-to-end benchmark

Because thinning gives a known set of newly created zeros, calculate:

- sensitivity/recall;
- precision/positive predictive value;
- F1 score;
- false-positive rate among zeros that were already zero before thinning;
- event q-value/confidence calibration where applicable.

Only after detection is scored should expression-recovery metrics be calculated.

Report at least two recovery quantities:

1. recovery error among all true artificial technical zeros;
2. recovery error conditional on the event being correctly detected.

The second isolates the magnitude estimator; the first reflects the practical end-to-end pipeline.

---

# 21. Reference value for count-thinning recovery

If recovery operates on normalized/log expression, the reference magnitude should be the value obtained by applying the same normalization pipeline to the **unthinned original counts**.

Define

\[
X^{ref}=Normalize(Y)
\]

and

\[
X^{thin}=Normalize(Y').
\]

For a thinning-created zero, compare the recovered value from the thinned pipeline with \(X^{ref}_{gc}\).

Do not compare a normalized recovery value directly with raw counts.

---

# 22. Paired statistical comparison

Artificial events within the same gene, membership, and cell state are correlated. Therefore a naive standard error treating every coordinate as independent is too optimistic.

Use paired resampling at a higher level, preferably:

- benchmark seed;
- gene;
- biological stratum;

or a hierarchical bootstrap.

For a metric \(L\), compare two methods by

\[
\Delta L_s
=L_{A,s}-L_{B,s}
\]

on the same seed/mask.

Report:

- median paired difference;
- mean paired difference;
- bootstrap 95% confidence interval;
- fraction of seeds in which each method wins.

---

# 23. Recommended primary decision criteria

`tree_local_factor` should not become the default merely because one aggregate metric improves.

A defensible promotion rule is:

1. lower paired RMSE than `masked_factor` on the locked oracle benchmark;
2. lower paired RMSE than the tree+embedding local-mean ablation, demonstrating incremental factor value;
3. no meaningful increase in absolute bias;
4. higher or comparable concordance correlation;
5. better or comparable gene-level variance/distribution preservation;
6. 95% predictive coverage reasonably close to nominal without extreme interval inflation;
7. no systematic degradation in any major trusted biological state;
8. the advantage persists under count-level thinning.

If tree-local weighting improves only the oracle benchmark but fails after thinning, it should remain an optional recovery engine rather than the package default.

---

# 24. Minimal implementation sketch

The benchmark should be implemented independently of the recovery internals so that benchmark code cannot accidentally share hidden state with the estimator.

A minimal oracle experiment is conceptually:

```r
set.seed(seed)

X_ref <- x
eligible <- which(X_ref > 0, arr.ind = TRUE)

# sample a constrained subset of positive coordinates
chosen <- sample(seq_len(nrow(eligible)), size = floor(0.10 * nrow(eligible)))
idx <- eligible[chosen, , drop = FALSE]

X_mask <- X_ref
X_mask[idx] <- 0

mask <- Matrix::sparseMatrix(
  i = idx[, 1],
  j = idx[, 2],
  x = TRUE,
  dims = dim(X_ref),
  dimnames = dimnames(X_ref)
)

truth <- X_ref[idx]

fit <- recover_dropout_expression(
  x = X_mask,
  mask = mask,
  membership = membership_fit,
  embedding = embedding_ref,
  recovery_method = "tree_local_factor",
  return_details = TRUE
)

pred <- fit$events$recovered
rmse <- sqrt(mean((pred - truth)^2))
mae <- mean(abs(pred - truth))
bias <- mean(pred - truth)
```

Production benchmark code must add the stratified masking constraints described above rather than relying on the simple unconstrained `sample()` shown here.

---

# 25. Required benchmark outputs

Each run should write an event-level table containing at least:

```text
seed
mask_fraction
method
gene
cell
hard_stratum
membership
truth
prediction
error
absolute_error
original_expression_bin
gene_detection_rate
positive_donor_count
effective_donors
tree_distance_weighted_mean
embedding_distance_weighted_mean
predictability
shrinkage
prediction_sd
recovery_method
```

This event table is the audit trail from which all summary figures and statistics should be reproducible.

---

# 26. Required summary tables

At minimum generate:

## Overall

```text
method | RMSE | MAE | Bias | Pearson | Spearman | CCC
```

## By expression level

```text
method | expression_bin | n | RMSE | MAE | Bias
```

## By effective donor support

```text
method | n_eff_bin | n | RMSE | MAE | prediction_sd | coverage_95
```

## By biological state

```text
method | state | n | RMSE | MAE | Bias | CCC | variance_error
```

## By tree distance

```text
method | tree_distance_bin | n | RMSE | MAE | prediction_sd | coverage_95
```

---

# 27. Required figures

Recommended figures:

1. truth vs recovered scatter with identity line;
2. per-method RMSE/MAE paired across seeds;
3. error versus original expression magnitude;
4. error versus effective donor number;
5. error versus tree distance;
6. prediction SD versus absolute error;
7. nominal versus empirical interval coverage;
8. original versus recovered gene variance;
9. Wasserstein distance by method;
10. biological-state-specific RMSE heatmap.

The truth-vs-recovered plot must include the identity line

\[
y=x
\]

because correlation alone can conceal systematic shrinkage or inflation.

---

# 28. Benchmark invariants

The benchmark itself must enforce the following invariants.

## 28.1 Hidden truth is never available to the estimator

The reference value must be stored outside `X_mask` before recovery.

## 28.2 All methods use exactly the same masks

No method-specific resampling.

## 28.3 Observed non-masked coordinates must remain unchanged

For every coordinate with artificial mask equal to zero:

\[
X^{out}_{gc}=X^{mask}_{gc}.
\]

Any violation is an automatic benchmark failure.

## 28.4 Target leakage remains prohibited

Artificially masked target genes must be treated exactly as real recovery targets. They may not re-enter the factor feature set simply because their original values are known to the benchmark harness.

## 28.5 Hyperparameter selection is separated from final testing

Locked test masks are never used to choose the recovery configuration.

---

# 29. Interpretation

The primary random-knockout experiment estimates

\[
\boxed{
E\left[
L(\widehat X_{gc},X^{ref}_{gc})
\mid
X^{ref}_{gc}>0
\right]
}
\]

for known missing positive expression.

It therefore directly tests the stated recovery target:

\[
E[X_{gc}\mid X_{gc}>0,\text{biological state}].
\]

It does **not** by itself prove that naturally observed zeros are technical dropouts. That is a detection question and belongs to the thinning/end-to-end benchmark.

This distinction is essential:

```text
oracle positive knockout
    -> validates expression magnitude recovery

count thinning + de novo detection
    -> validates technical-zero detection and complete workflow
```

---

# 30. Recommended execution order

Run the benchmark in this order:

```text
Phase A
random positive-coordinate knockout
fixed embedding + fixed hierarchy
oracle mask
        |
        v
component ablation
        |
        v
select candidate hyperparameters on development masks

Phase B
locked random-positive confirmatory masks
        |
        v
final recovery comparison
        |
        v
DV + uncertainty calibration

Phase C
raw-count binomial thinning
recompute complete state representation
        |
        v
detection + recovery benchmark
        |
        v
biological-state stress tests
```

The package default should be changed only after Phase B and Phase C support the same conclusion.
