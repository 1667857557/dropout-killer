# DropoutKiller 0.4.0

- Replaced neighbor-only recovery as the default with membership-local masked coexpression-factor recovery.
- Dropout-mask entries are treated as missing during recovery learning; unmasked zeros remain observed biological zeros.
- Excluded all current recovery-target genes from membership factor features so target expression cannot leak back into the cell-state representation used to predict it.
- Added target-gene ridge prediction with exact analytic leave-one-out shrinkage toward the membership mean; unsupported coexpression automatically collapses to the mean without explicit fold CV.
- Added event-level predictive standard deviations and a sparse `predictive_variance` output so uncertainty is not discarded when recovered means are written to the expression matrix.
- Added `sample_dropout_expression()` for uncertainty-aware completed-matrix draws; observed values remain exact.
- Kept `weighted_neighbor_prediction()` and `recovery_method = "neighbor"` as an explicit comparator/legacy engine, but it is no longer the default recovery path.
- Added synthetic masked-expression tests for coexpression recovery, fallback behavior, selective invariants, and uncertainty-aware draws.

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
