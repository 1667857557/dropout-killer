# Dropout detector source-parity contract

This note records the exact upstream mathematical contracts implemented by the
high-level DropoutKiller detector options in version 0.8.

## `scgacl_gamma_normal` (default)

Upstream implementation: `Hhjyl/scGACL`,
`evaluation/dropout_identify.py` and `trainer/config.py`.

The released configuration uses `IDENTIFY_USE_CELLTYPE = True` and
`DROP_THRESHOLD = 0.5`. DropoutKiller therefore requires explicit `group`
subpopulation labels for this source-faithful path rather than silently replacing
the upstream Scanpy/sklearn fallback with a different R clustering implementation.

For raw count `C`, scGACL first defines

```text
x = log(1.01 + C)
point = log(1.01)
```

and fits one mixture per gene and supplied cell subpopulation:

```text
f(x) = lambda Gamma(x; alpha, beta)
       + (1-lambda) Normal(x; mu, sigma)
```

where `beta` is a Gamma rate. The Gamma responsibility is

```text
d(x) = lambda Gamma(x; alpha, beta) / f(x).
```

Initialization, invalid-gene rules, EM updates, population standard deviation,
Gamma-shape equation, base-10 log-likelihood stopping rule, and the released
iteration limit are reproduced in `R/scgacl_detection.R`. Only raw-zero
coordinates are eligible for the DropoutKiller mask. The released `main.py`
retains reconstructed values unless `predict_drop < DROP_THRESHOLD`; therefore
the code-faithful boundary is `d >= 0.5`.

For a fixed gene and cell subpopulation every raw zero has the same transformed
value `point`, so evaluating `d(point)` once and assigning it to all zero
coordinates is exactly equivalent to building the upstream dense posterior matrix.
The R implementation processes genes in bounded batches; batching and sparse
storage are tested to leave event coordinates and posterior values unchanged.

A frozen cross-language regression test reproduces the released Python/SciPy
parameters and zero posterior for a fixed count vector.

## `alra_global` (optional comparator)

Upstream implementation: `KlugerLab/ALRA`, `R/alra.R`.

DropoutKiller uses one all-cell matrix. No `group`, `split_by`, or SuperCell
membership fragments the ALRA SVD.

For cells by genes matrix `A`, the upstream normalization is

```text
A_cg <- log(1 + 10000 * C_cg / sum_h C_ch).
```

Automatic rank selection reproduces:

```text
K = 100
noise_start = 80
q_choose = 2
threshold = 6 SD
```

with the original source validity checks (`K <= min(dim(A))` and at least five
noise singular values). The final randomized SVD uses `q = 10`.

The dropout-detection event set is the zero-coordinate support of ALRA's adaptive
thresholding stage:

```text
tau_g = abs(Q_0.001(Ahat_.g))
event_cg = (A_cg == 0) & (Ahat_cg > tau_g).
```

This is the exact zero-identification gate in the original `alra()` function.
DropoutKiller intentionally does not use ALRA's later mean/SD rescaling to recover
expression magnitudes, because recovery is a separate P1 stage. A source-parity
unit test recomputes the upstream randomized SVD and quantile gate directly and
requires identical event coordinates, thresholds, and low-rank values.

The historical option name `alra_global_by_group` is deprecated and maps to
`alra_global`; it no longer fits cell-class-specific SVDs.
