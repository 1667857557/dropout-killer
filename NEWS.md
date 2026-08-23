# DropoutKiller 0.2.0

- Rebuilt membership construction around the SuperCell graph/coarse-graining contract:
  biological hard strata -> embedding-space kNN graph -> walktrap hierarchy -> gamma-controlled memberships.
- Added membership-local ALRA-derived zero gating and negative-tail confidence scoring.
- Added event-only Gaussian neighbor borrowing and sparse selective replacement.
- Restricted biological priors to **PPI and pathway information only**. GRN/TF-target priors are intentionally unsupported.
- Added membership-local standardized PPI/pathway prediction and explicit PPI/pathway fusion weights.
- Added Seurat integration, diagnostics, mathematical contract, unit tests, and GitHub Actions R CMD check.
- Enforced zero-only masks and exact preservation of observed non-zero expression.
