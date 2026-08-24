# DropoutKiller 0.3.0

- Removed the PPI/pathway prior branch completely, including prior constructors, preparation/fusion/prediction APIs, internal prior prediction code, documentation, and tests.
- Removed `alpha`, `ppi`, `pathway`, `prior_weights`, `gene_prior`, and `prior_standardize` from recovery APIs.
- Recovery now uses only membership-constrained Gaussian neighbor borrowing; the effective donor weights are estimated from latent-space distances and renormalized over eligible donors.
- Removed prior-specific event fields and settings while preserving zero-only replacement, sparse output, and exact preservation of observed non-zero values.

# DropoutKiller 0.2.0

- Rebuilt membership construction around the SuperCell graph/coarse-graining contract:
  biological hard strata -> embedding-space kNN graph -> walktrap hierarchy -> gamma-controlled memberships.
- Added membership-local ALRA-derived zero gating and negative-tail confidence scoring.
- Added event-only Gaussian neighbor borrowing and sparse selective replacement.
- Added Seurat integration, diagnostics, mathematical contract, unit tests, and GitHub Actions R CMD check.
- Enforced zero-only masks and exact preservation of observed non-zero expression.
