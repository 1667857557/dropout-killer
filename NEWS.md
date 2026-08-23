# DropoutKiller 0.2.0

- Replaced the incorrect embedding-to-`SCimplify()` wrapper with an internal SuperCell-style kNN/walktrap membership engine.
- Added pre-graph biological strata, approximate anchor/centroid assignment, and membership diagnostics.
- Replaced the arbitrary sigmoid score with ALRA-gated negative-tail confidence scoring restricted to observed zeros.
- Added sparse mask construction and membership-local SVD processing.
- Added Gaussian latent-space neighbor borrowing using positive-expression donors.
- Added scale-aware PPI/pathway/GRN prediction and multi-network combination.
- Added selective 75/25 hybrid recovery with availability-aware weight renormalization.
- Added Seurat integration, validation invariants, documentation, unit tests, and R CMD check CI.
