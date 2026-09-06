# DropoutKiller 0.8 mathematical contract

## 1. Detection and recovery are different statistical problems

Let raw gene-by-cell counts be

\[
C=(C_{gc})\in\mathbb R_+^{G\times N}.
\]

DropoutKiller never treats every observed zero as technical dropout. It first constructs a sparse detector mask \(M\), then estimates positive expression only at coordinates selected by \(M\). Observed non-selected coordinates are immutable on the recovery working scale.

Version 0.8 deliberately separates the detector input from the recovery input:

```text
raw counts C
  |-- scGACL detector scale: log(1.01 + C)
  |
  `-- recovery scale: library normalize to 1e4 -> log1p
```

This prevents a detector-specific transformation from silently changing the validated P1 recovery model.

---

## 2. Production detector: released scGACL Gamma-Normal mixture

The default is

```r
detection_method = "scgacl_gamma_normal"
scgacl_dropout_threshold = 0.5
```

The implementation is a direct mathematical port of
`Hhjyl/scGACL/evaluation/dropout_identify.py`.

### 2.1 Transform

The released detector transforms raw counts as

\[
x_{gc}=\log(1.01+C_{gc}),
\]

so every raw zero is mapped to

\[
x_0=\log(1.01).
\]

No ALRA library normalization is applied before scGACL detection.

### 2.2 Subpopulation scope

The released scGACL configuration sets `IDENTIFY_USE_CELLTYPE = True`. The
source-faithful DropoutKiller implementation therefore requires `group`, normally
major cell type, and fits each gene independently inside each supplied group.

The package does not silently replace scGACL's alternative Python/Scanpy
preclustering route with a different R clustering implementation.

### 2.3 Mixture likelihood

For one gene inside one subpopulation:

\[
f(x)=\lambda f_\Gamma(x;\alpha,\beta)
+(1-\lambda)f_N(x;\mu,\sigma),
\]

where \(\beta\) is the Gamma **rate**.

The Gamma responsibility is

\[
w_i=
\frac{\lambda f_\Gamma(x_i;\alpha,\beta)}
{\lambda f_\Gamma(x_i;\alpha,\beta)+(1-\lambda)f_N(x_i;\mu,\sigma)}.
\]

For an observed zero, the scGACL dropout score is simply this responsibility
evaluated at \(x_0\):

\[
d_{g,t}=w(x_0).
\]

Because every raw zero in the same gene/subpopulation has exactly \(x=x_0\), all
such zero coordinates receive the same posterior. Evaluating \(d_{g,t}\) once is
algebraically identical to generating the released full dense posterior matrix.

### 2.4 Released initialization and gene guards

For a gene vector \(x\):

\[
\lambda^{(0)}=\frac{\#\{i:x_i=x_0\}}{n}.
\]

The released code:

- marks a fit invalid when \(\lambda^{(0)}>0.95\);
- replaces \(\lambda^{(0)}=0\) by 0.01;
- initializes \(\alpha=1.5\), \(\beta=1\);
- initializes \(\mu\) and the **population** standard deviation \(\sigma\) from values \(x>x_0\);
- replaces an initial zero standard deviation by 0.01;
- marks a gene invalid before EM when
  \[
  |\bar x-x_0|<10^{-2}.
  \]

These details are retained because changing them changes the detector.

### 2.5 EM updates

With Gamma responsibility \(w_i\):

\[
\lambda\leftarrow\frac{1}{n}\sum_i w_i,
\]

\[
\mu\leftarrow
\frac{\sum_i(1-w_i)x_i}{\sum_i(1-w_i)},
\]

\[
\sigma\leftarrow
\sqrt{
\frac{\sum_i(1-w_i)(x_i-\mu)^2}{\sum_i(1-w_i)}
}.
\]

For the Gamma component define

\[
S=\sum_iw_i,
\qquad
T=\sum_iw_ix_i,
\qquad
U=\sum_iw_i\log x_i,
\]

and

\[
v=-\frac{U}{S}-\log\frac{S}{T}.
\]

If \(v\le0\), the released code sets \(\alpha=20\). Otherwise it initializes

\[
\alpha_0=
\frac{3-v+\sqrt{(v-3)^2+24v}}{12v},
\]

caps \(\alpha\) at 20 when \(\alpha_0\ge20\), or solves the unique positive root

\[
\log\alpha-\psi(\alpha)=v.
\]

Then

\[
\beta=\frac{S}{T}\alpha.
\]

The released stopping criterion is the squared change in base-10 log likelihood:

\[
\left(\ell_{10}^{(r)}-\ell_{10}^{(r-1)}\right)^2\le0.5,
\]

or the released iteration limit.

### 2.6 Zero call

The paper writes the decision as a posterior threshold. The released `main.py`
retains a zero unless `predict_drop < threshold`; therefore the exact code-boundary
behavior is

\[
M_{gc}=1
\iff
C_{gc}=0
\ \text{and}\ 
d_{g,t(c)}\ge0.5.
\]

---

## 3. Optional original ALRA detector

The source-faithful comparator is

```r
detection_method = "alra_global"
```

It is **one all-cell ALRA**, not one ALRA per cell class.

### 3.1 Original ALRA normalization

For library size

\[
L_c=\sum_g C_{gc},
\]

ALRA uses

\[
A_{cg}=\log\left(1+10^4\frac{C_{gc}}{L_c}\right),
\]

with cells in rows and genes in columns.

### 3.2 Original automatic rank

The released implementation computes up to \(K=100\) singular values with a
randomized SVD using \(q=2\) power iterations. Consecutive spacings are

\[
\Delta_i=\sigma_i-\sigma_{i+1}.
\]

With `noise_start = 80`, the source estimates the noise-spacing mean and standard
deviation from `diffs[noise_svals - 1]`, then selects

\[
k=\max\left\{i:
\frac{\Delta_i-\mu_\Delta}{s_\Delta}>6
\right\}.
\]

The original validity constraints are retained:

- `K` may not exceed the smaller matrix dimension;
- at least five singular values must be assigned to the noise region.

DropoutKiller no longer modifies these rules to force ALRA to run on small
cell-class blocks.

### 3.3 Original low-rank gate

The final rank-\(k\) randomized SVD uses \(q=10\):

\[
\widehat A=U_kD_kV_k^T.
\]

For every gene:

\[
\tau_g=
\left|Q_{0.001}(\widehat A_{\cdot g})\right|,
\]

using R's default empirical quantile rule as in the released code. An observed zero
is selected when

\[
A_{cg}=0
\quad\text{and}\quad
\widehat A_{cg}>\tau_g.
\]

The old public value `alra_global_by_group` is deprecated and maps to this all-cell
implementation with a warning.

---

## 4. Historical local detectors

For reproducibility only:

```r
detection_method = "eb_zero_null"
detection_method = "alra_quantile"
```

These retain the pre-0.8 membership-local detector implementations and are not the
production default.

---

## 5. Membership construction and recovery geometry

The final SuperCell membership is a recovery borrowing block, not an scGACL or
ALRA estimation block. Supplied `group` and `split_by` may further impose hard
recovery strata. Within a final membership, the retained walktrap hierarchy and
biological embedding provide continuous donor weights.

The default recovery remains

```r
recovery_method = "p1_stabilized_state"
```

For a target-gene fold, the fold is excluded from predictor construction. The
standardized non-target predictor state receives one row-stochastic geometry
smoothing step

\[
Z_{stable}=(1-\rho)Z+\rho PZ,
\qquad \rho=0.25,
\]

followed by target-safe SVD and positive-donor ridge recovery. Analytic leave-one-
out diagnostics determine shrinkage and bias calibration.

Observed non-selected coordinates are not changed.

---

## 6. Provenance and verification contract

### scGACL

Published implementation:

<https://github.com/Hhjyl/scGACL>

Detector source:

<https://github.com/Hhjyl/scGACL/blob/main/evaluation/dropout_identify.py>

DropoutKiller regression tests freeze mixture parameters and zero posterior values
computed by that released Python implementation. The R port must reproduce them
within numerical tolerance.

### ALRA

Original implementation:

<https://github.com/KlugerLab/ALRA>

The DropoutKiller comparator follows the released `normalize_data()`, `choose_k()`
and `alra()` mathematical path. Tests explicitly reject the former adaptive
small-cell-class behavior.
