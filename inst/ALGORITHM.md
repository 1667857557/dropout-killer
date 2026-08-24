# DropoutKiller 0.5 mathematical contract

## 1. Problem definition

For one biological membership, let

\[
X=(X_{gc})\in\mathbb R_+^{G\times n}
\]

be the supplied normalized expression matrix. DropoutKiller does not claim that a zero itself identifies technical dropout. It separates the problem into:

1. **detection**: decide which observed zeros have sufficient evidence against a biological-zero reconstruction null;
2. **recovery**: conditional on an event already being classified as technical dropout, estimate its positive latent expression magnitude and uncertainty.

Observed non-dropout coordinates are immutable.

If the final mask is \(M\), then

\[
X^{out}_{gc}=
\begin{cases}
X_{gc}, & M_{gc}=0,\\
\widehat\lambda_{gc}, & M_{gc}=1.
\end{cases}
\]

The package never writes a prediction over an observed nonzero value.

---

## 2. Membership construction

Cells may first be partitioned by supplied hard biological strata, for example major cell type and optionally donor/condition. Within each stratum, a Euclidean kNN graph is built from the supplied low-dimensional embedding and coarse-grained with a SuperCell-style target size controlled by `gamma`.

Default:

\[
\gamma=150.
\]

Membership labels are canonicalized once in cell order before returning. This changes labels only, not the underlying partition, and makes later alignment idempotent.

---

## 3. Membership-local low-rank reconstruction

For membership \(m\), let

\[
X_m\in\mathbb R_+^{G\times n_m}.
\]

A rank-\(k\) reconstruction is obtained from a truncated SVD:

\[
X_m\approx U_kD_kV_k^T=\widehat X_m.
\]

Automatic rank selection keeps the existing singular-spacing heuristic. The 0.5 redesign changes the **zero-null thresholding**, not the low-rank reconstruction itself, so detection and recovery changes can be attributed separately.

---

## 4. Why the original empirical 0.1% threshold is not the default

The historical ALRA-style gate used a per-gene empirical lower quantile

\[
Q_{0.001}(\widehat X_{g\cdot}).
\]

In a membership of \(n_m\) cells, the expected number of observations represented by this tail is

\[
0.001n_m.
\]

For \(n_m<1000\), the empirical 0.1% quantile is determined almost entirely by the first one or two order statistics. This creates an avoidable finite-sample instability in membership-local use.

The historical implementation is retained under

```r
detection_method = "alra_quantile"
```

but is not the 0.5 default.

---

## 5. Finite-sample empirical-Bayes biological-zero null

The default detector is

```r
detection_method = "eb_zero_null"
```

and retains the ALRA symmetry idea:

> under a biological-zero null, low-rank reconstruction error is approximately symmetric around zero.

For gene \(g\), define negative reconstructed values

\[
\mathcal N_g=\{c:\widehat X_{gc}<0\}
\]

with size

\[
n_g^-=|\mathcal N_g|.
\]

The local second-moment estimator is

\[
s_g^2=
\frac{1}{\max(n_g^-,1)}
\sum_{c\in\mathcal N_g}\widehat X_{gc}^2.
\]

Because this estimate is noisy when \(n_g^-\) is small, define a robust membership-level prior center

\[
s_0^2=\operatorname{median}_{g\in\mathcal G_*}(s_g^2),
\]

where \(\mathcal G_*\) contains genes with at least `min_negative` negative reconstructed values and finite positive variance.

Let

\[
\nu_0=\texttt{variance_prior_df}.
\]

The shrinkage weight is

\[
w_g=\frac{n_g^-}{n_g^-+\nu_0},
\]

and the finite-sample zero-null variance is

\[
\boxed{
\widetilde s_g^2=w_gs_g^2+(1-w_g)s_0^2
}
\]

with

\[
\widetilde s_g=\sqrt{\widetilde s_g^2}.
\]

Thus:

- large memberships / many negative residuals: \(w_g\to1\), gene-specific evidence dominates;
- small memberships / weak negative support: \(w_g\to0\), the estimate shrinks toward the membership-level scale;
- finite-sample uncertainty therefore grows smoothly rather than relying on an extreme empirical order statistic.

---

## 6. Zero-event hypothesis testing

Only observed zeros with positive low-rank reconstruction are candidates for the expressed/dropout alternative.

For candidate \((g,c)\):

\[
Z_{gc}=\frac{\widehat X_{gc}}{\widetilde s_g}.
\]

Under the working one-sided biological-zero null,

\[
p_{gc}=P(N(0,1)\ge Z_{gc})
       =1-\Phi(Z_{gc}).
\]

Within each **gene and membership**, candidate p-values are Benjamini-Hochberg adjusted while using the full number of observed zeros of that gene as the number of tested hypotheses. Let the adjusted value be

\[
q_{gc}.
\]

DropoutKiller stores

\[
\boxed{
\text{confidence}_{gc}=1-q_{gc}
}
\]

for compatibility with the existing score interface.

This is **not** a Bayesian posterior probability of technical dropout.

With the default

```r
threshold = 0.95
```

selection corresponds to

\[
q_{gc}\le0.05
\]

under the working zero-null model.

The finite-sample detector therefore converts the old score cutoff into an interpretable gene-wise multiple-testing threshold rather than an uncalibrated probability label.

---

## 7. Detection and recovery are separate estimands

A key 0.5 distinction is:

\[
\text{detect whether a zero is incompatible with biological zero}
\]

is different from

\[
\text{estimate its expression magnitude once classified as dropout}.
\]

The recovery model never changes the mask produced by detection.

---

## 8. Target-leakage-free membership factor state

For a membership containing recovery events, let \(T\) be the set of all target genes with at least one masked event.

Every gene in \(T\) is excluded from factor-feature learning.

Using a set \(F\subseteq\{1,\ldots,G\}\setminus T\) of high-variance non-target genes, standardized expression is decomposed by truncated SVD to obtain cell factor scores

\[
z_c\in\mathbb R^K.
\]

Therefore the prediction path for a target \((g,c)\) is

\[
X_{c,-T}\rightarrow z_c\rightarrow\widehat X_{gc}
\]

and never

\[
X_{gc}\rightarrow z_c\rightarrow\widehat X_{gc}.
\]

This prevents direct target leakage.

---

## 9. Positive-conditional target magnitude

Once \((g,c)\) is already classified as technical dropout, the default recovery estimand is

\[
\boxed{
E[X_{gc}\mid X_{gc}>0,z_c,m]
}
\]

rather than

\[
E[X_{gc}\mid z_c,m]
\]

that mixes reliable biological zeros with positive expression magnitudes.

This is controlled by

```r
factor_target = "positive"
```

(default).

For reproducibility, the previous unconditional target remains available as

```r
factor_target = "all_observed"
```

Unmasked zeros are never changed; they are simply excluded from the **magnitude model** after a separate detector has already classified the query coordinate as technical dropout.

This design deliberately couples a more selective FDR-aware detector to a positive-conditional recovery model. Using positive-conditional recovery on arbitrary natural zeros would be biased upward and is not supported.

---

## 10. Ridge factor regression

For target gene \(g\), let \(D_g\) be reliable donor cells. Under the default positive target,

\[
D_g=\{c:M_{gc}=0,\ X_{gc}>0\}.
\]

Construct

\[
Z_g=[\mathbf 1,z_c]_{c\in D_g}.
\]

The ridge estimator is

\[
\widehat\beta_g=
\arg\min_\beta
\left[
\|y_g-Z_g\beta\|_2^2+
\lambda\|\beta_{factor}\|_2^2
\right],
\]

with the intercept unpenalized.

Closed form:

\[
\widehat\beta_g=
(Z_g^TZ_g+P)^{-1}Z_g^Ty_g.
\]

---

## 11. Exact analytic leave-one-out shrinkage

For ridge linear smoother \(H\), the leave-one-out prediction is

\[
\widehat y_i^{(-i)}
=
y_i-
\frac{y_i-\widehat y_i}{1-h_{ii}}.
\]

Let the leave-one-out null prediction be the donor mean excluding cell \(i\):

\[
\mu_i^{(-i)}.
\]

Define

\[
d_i=\widehat y_i^{(-i)}-\mu_i^{(-i)},
\]

\[
t_i=y_i-\mu_i^{(-i)}.
\]

The squared-error-optimal linear shrinkage of the factor contribution is

\[
\boxed{
q_g=
\operatorname{clip}_{[0,1]}
\frac{\sum_i d_it_i}{\sum_i d_i^2}
}
\]

and query prediction is

\[
\boxed{
\widehat\lambda_{gc}
=
\max\{0,\mu_g+q_g(\widehat y^{factor}_{gc}-\mu_g)\}
}
\]

where \(\mu_g\) is the positive-donor mean under the default target.

If held-out factor information is unsupported,

\[
q_g=0
\]

and recovery falls back to the positive membership mean.

---

## 12. Predictability diagnostic

Define leave-one-out SSE for the donor-mean null and shrunken factor model:

\[
SSE_0=\sum_i(y_i-\mu_i^{(-i)})^2,
\]

\[
SSE_1=\sum_i(y_i-\widehat y_{i,shrunk}^{(-i)})^2.
\]

Stored predictability is

\[
\boxed{
D_g^2=
\operatorname{clip}_{[0,1]}
\left(1-\frac{SSE_1}{SSE_0}\right)
}
\]

when the null SSE is positive.

This measures whether cell-specific coexpression adds held-out information beyond the membership-level positive mean.

---

## 13. Predictive uncertainty and differential variability

The previous implementation estimated target residual variance from in-sample fitted residuals. Version 0.5 instead uses leave-one-out residual MSE:

\[
\widehat\sigma_{g,LOO}^2
=
\frac1{|D_g|}
\sum_{i\in D_g}
(y_i-\widehat y_{i,shrunk}^{(-i)})^2.
\]

For query cell \(c\), parameter uncertainty is approximated using ridge leverage \(h_c\) and mean-estimation contribution:

\[
\ell_c=
\frac{(1-q_g)^2}{|D_g|}
+q_g^2h_c.
\]

The stored predictive variance is

\[
\boxed{
V_{gc}=
\widehat\sigma_{g,LOO}^2(1+\ell_c)
}
\]

and

\[
\text{prediction\_sd}_{gc}=\sqrt{V_{gc}}.
\]

This is an approximate predictive variance, not a fully calibrated Bayesian posterior.

For the default positive target, repeated completed draws use a Gamma moment match. Given stored predictive mean \(m_{gc}>0\) and variance \(V_{gc}>0\), define

\[
k_{gc}=\frac{m_{gc}^2}{V_{gc}},
\qquad
\theta_{gc}=\frac{V_{gc}}{m_{gc}}.
\]

Then

\[
X_{gc}^{(b)}\sim\operatorname{Gamma}(k_{gc},\theta_{gc})
\]

satisfies

\[
E[X_{gc}^{(b)}]=m_{gc},
\qquad
\operatorname{Var}(X_{gc}^{(b)})=V_{gc},
\qquad
X_{gc}^{(b)}>0.
\]

This preserves the first two predictive moments exactly under the approximation and avoids the mean shift produced by truncating a Gaussian below zero.

---

## 14. Why the mean matrix alone cannot preserve DV

For latent expression \(\lambda_{gc}\):

\[
\boxed{
\operatorname{Var}(\lambda_g\mid Y)
=
\operatorname{Var}_c(E[\lambda_{gc}\mid Y])
+
E_c(\operatorname{Var}[\lambda_{gc}\mid Y])
}
\]

Replacing a missing value by only its conditional mean discards the second term and therefore contracts variance.

Accordingly, DropoutKiller exposes:

- deterministic recovery mean;
- event-level `prediction_sd`;
- sparse `predictive_variance`;
- `sample_dropout_expression()` for repeated completed draws.

DV/covariance/network analysis should propagate repeated draws rather than treating the deterministic mean matrix as error-free observation.

---

## 15. Multiple-imputation interpretation

For dropout events \(\mathcal D\), generate

\[
X^{(1)},\ldots,X^{(B)}.
\]

Observed coordinates are identical across all draws. Only \((g,c)\in\mathcal D\) vary.

For a statistic \(T\), uncertainty-aware downstream analysis should operate on

\[
T^{(b)}=T(X^{(b)})
\]

and summarize the distribution over \(b\), rather than computing only

\[
T(E[X\mid Y]).
\]

---

## 16. Validation contract

Recovery cannot validate itself by showing that post-recovery coexpression becomes stronger. If the predictor used coexpression to generate \(\widehat X_g\), an increase in

\[
\operatorname{cor}(\widehat X_g,X_h)
\]

may merely reflect model feedback.

Validation must use held-out information.

### 16.1 Pseudo-mask validation

Hide reliable observed positive coordinates, fit without them, and compare held-out predictions against their known observations.

### 16.2 Count-level thinning

For UMI counts, prefer

\[
Y'_{gc}\mid Y_{gc}\sim\operatorname{Binomial}(Y_{gc},\rho),
\]

or equivalent Poisson thinning. Newly created zeros then have known technical origin and better approximate the measurement process than setting normalized positives directly to zero.

### 16.3 Separate detection and recovery metrics

Report separately:

- detection recall / mask stability;
- recovery error conditional on an oracle mask;
- end-to-end error;
- predictive interval coverage;
- DV/covariance recovery.

---

## 17. Statistical boundaries

1. `confidence` under `eb_zero_null` is `1 - BH q`, not posterior dropout probability.
2. Gene-wise BH control is conditional on the approximate symmetric Gaussian zero-null and does not constitute exact global FDR control across all genes and memberships.
3. Positive-conditional recovery is appropriate only for coordinates already selected as technical dropout; it must not be applied indiscriminately to natural zeros.
4. Factor rank describes only the predictable coexpression component. Residual biological variability remains in predictive variance.
5. The historical empirical-quantile detector and neighbor recovery engine remain available only as explicit reproducibility/comparator paths.
6. Recovered values are continuous normalized-expression estimates, not integer raw counts.

---

## 18. Default 0.5 workflow

```text
hard biological strata
        |
        v
SuperCell-style membership
        |
        v
membership low-rank reconstruction
        |
        v
negative-null variance EB shrinkage
        |
        v
one-sided p values + gene-wise BH
        |
        v
q <= 0.05 dropout mask
        |
        v
exclude target genes from factor features
        |
        v
positive-donor target ridge
        |
        v
analytic LOO shrinkage
        |
        v
mean + predictive variance
        |
        v
selective replacement + optional repeated draws
```
