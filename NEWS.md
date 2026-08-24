# DropoutKiller 0.5.0

- Replaced the default membership-local empirical `quantile_prob = 0.001` ALRA gate with a finite-sample empirical-Bayes zero-null detector. The historical quantile gate remains available with `detection_method = "alra_quantile"`.
- The new default `detection_method = "eb_zero_null"` estimates gene-specific zero-null variance from negative low-rank reconstructions and shrinks it toward a robust membership-level variance center using `variance_prior_df` pseudo-observations.
- Positive reconstructed values at observed zeros are tested against the symmetric zero-null with one-sided Gaussian p-values and gene-wise Benjamini-Hochberg correction. Event `confidence` is now `1 - q_value`, not a Bayesian posterior probability; the default `threshold = 0.95` therefore corresponds to gene-wise q <= 0.05 under the working null model.
- Changed masked-factor recovery to estimate positive conditional magnitude after an event has already been classified as technical dropout. The default `factor_target = "positive"` fits target genes only on reliable positive donors; `factor_target = "all_observed"` reproduces the previous zero-inflated target.
- Replaced in-sample residual variance in the factor engine with analytic leave-one-out residual MSE. Positive-conditional multiple-imputation draws now use Gamma moment matching so the sampled distribution remains non-negative while preserving the stored predictive mean and variance.
- Canonicalized internally built membership labels before returning them so returned membership vectors, detection labels, and recovery event labels remain idempotently aligned without changing the partition.
- Preserved historical public argument positions by appending `detection_method`, `variance_prior_df`, and `factor_target` after the existing 0.4 API.
- Added regression tests for finite-sample variance shrinkage, gene-wise BH score semantics, positive-conditional recovery, positive predictive draws, legacy quantile detection, and membership-label idempotence.

# DropoutKiller 0.4.0

- Replaced neighbor-only recovery as the default with membership-local masked coexpression-factor recovery.
- Dropout-mask entries are treated as missing target values during recovery; unmasked zeros remain observed biological zeros.
- Excluded all current recovery-target genes from membership factor features so target expression cannot leak back into the cell-state representation used to predict it.
- Simplified factor-state learning to one direct truncated SVD on standardized non-target genes; large memberships no longer build or fully diagonalize a dense cell-by-cell Gram matrix.
- Added target-gene ridge prediction with exact analytic leave-one-out shrinkage toward the membership mean; unsupported coexpression automatically collapses to the mean without explicit fold CV.
- Added event-level predictive standard deviations and a sparse `predictive_variance` output for masked-factor recovery so uncertainty is not discarded when recovered means are written to the expression matrix.
- Added `sample_dropout_expression()` for uncertainty-aware completed-matrix draws; observed values remain exact.
- Neighbor recovery now reports `uncertainty_available = FALSE` and `predictive_variance = NULL`; uncertainty-aware sampling rejects these results instead of interpreting missing uncertainty as zero variance.
- Preserved the pre-0.4 positional order of the legacy neighbor/recovery arguments and appended the new factor-engine options after the historical public API.
- Kept `weighted_neighbor_prediction()` and `recovery_method = "neighbor"` as an explicit comparator/legacy engine, but it is no longer the default recovery path.
- Added synthetic masked-expression, positional-compatibility, scaling-path, selective-invariant, and uncertainty-semantics tests.

# DropoutKiller 0.3.0

- Removed the PPI/pathway prior branch completely, including prior constructors, preparation/fusion/prediction APIs, internal prior prediction code, documentation, and tests.
- Removed `alpha`, `ppi`, `pathway`, `prior_weights`, `gene_prior`, and `prior_standardize` from recovery APIs.
- Recovery now uses only membership-constrained Gaussian neighbor borrowing; the effective donor weights are estimated from latent-space distances and renormalized over eligible donors.
- Removed prior-specific event fields and settings while preserving zero-only replacement, sparse output, and exact preservation of observed non-zero values.
- Changed the default SuperCell-style membership graining level from `gamma = 20` to `gamma = 150` in both membership construction and the end-to-end workflow.

# DropoutKiller 0.2.0

- Rebuilt membership construction around the SuperCell graph/coarse-graining contract:
  biological hard strata -> embedding-space kNN graph -> walktrap hierarchy -> gamma-controlled memberships.
- Added membership-local ALRA-derived zero gating and negative-tail confidence scoring.
- Added event-only Gaussian neighbor borrowing and sparse selective replacement.
- Added Seurat integration, diagnostics, mathematical contract, unit tests, and GitHub Actions R CMD check.
- Enforced zero-only masks and exact preservation of observed non-zero expression.
