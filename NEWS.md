# DropoutKiller 0.7.0

- Added `recovery_method = "barycentric"` as a new recovery-only candidate for the setting where dropout coordinates, biological embedding, SuperCell membership, and hierarchy are already available.
- Reinterpreted final membership as a soft prior for this engine rather than an absolute borrowing wall. Explicit hard biological strata remain strict zero-weight boundaries; embedding-local candidate cells from nearby sibling memberships may contribute.
- Added a hierarchy/embedding/membership geometry prior and a target-independent barycentric weight learner. For query cell `c`, donor weights solve a simplex-constrained convex reconstruction problem in biological-embedding space, regularized toward the geometry prior.
- Reused one learned cell-state weight vector across all dropout genes of the same query cell. Target-gene expression is excluded from weight learning; each target is recovered only from reliable positive donors after gene-specific weight renormalization.
- Added Kish effective positive-donor size, local positive prevalence, weighted biological variance, working predictive standard deviation, state-reconstruction improvement, same-membership weight, and tree/embedding distance diagnostics.
- Kept `masked_factor`, `tree_local_factor`, and `neighbor` unchanged as benchmark engines; `barycentric` is intentionally not promoted to the default before oracle-mask and thinning benchmarks establish a consistent advantage.
- Added a literature-and-design survey in `inst/LOCAL_RECOVERY_ARCHITECTURES.md` covering statistical smoothing, self-representation, Bayesian shrinkage, low-rank completion, graph methods, deep count models, local nonlinear regression, uncertainty propagation, and recovery-engine stacking.
- Added regression tests for simplex weights, soft cross-membership borrowing, absolute hard-stratum boundaries, selective zero-only writing, and uncertainty reporting.

# DropoutKiller 0.6.0

- Preserved the full `cluster_walktrap()` hierarchy inside `DropoutKillerMembership` instead of discarding it after the gamma cut. Exact builds retain per-stratum hierarchy and a cached tree index; approximate builds retain the anchor hierarchy and fall back to embedding-only weighting when a queried cell has no exact tree leaf.
- Added `recovery_method = "tree_local_factor"` directly to the original recovery dispatch. No public recovery function is overridden or redefined outside the existing call chain.
- Replaced the equal-weight positive membership fallback with a query-specific positive baseline. The final SuperCell membership is now the recovery borrowing block; walktrap-tree proximity and adaptive embedding distance determine continuous donor weights inside that block, while explicit hard strata may further split it.
- Added Kish effective donor size, weighted local variance, local positive prevalence, tree-distance and embedding-distance diagnostics for tree-local recovery events.
- Changed coexpression refinement to predict residual expression after subtracting leave-one-out local biological baselines. Donor influence is weighted by biological proximity to the actual dropout query cells.
- Reworked the initial event-wise tree-local ridge implementation for large high-confidence masks. Geometry now caches `W`, `W^2`, row weights, bandwidths, and cell-level distance summaries once; local gene statistics are evaluated in batches with matrix multiplication.
- Replaced one weighted ridge solve per dropout event with one residual ridge per target gene × final membership. Query-cell donor weights are aggregated for the gene-level fit, while `n_eff` keeps query-specific information shrinkage and uncertainty. Event count now affects indexing/output rather than the number of regression solves.
- Pre-indexed events by recovery block and target gene, removing repeated whole-event scans inside the gene loop.
- Added held-out factor shrinkage and effective-donor information shrinkage so weakly supported query neighborhoods collapse to the query-specific local mean rather than the whole-membership mean.
- Preserved the existing `masked_factor` and `neighbor` engines for component ablation and backward compatibility. Historical positional argument slots remain unchanged; tree-local controls are appended after the 0.5 API.
- Added regression tests for retained hierarchy, symmetric tree distance, membership-local borrowing, cached geometry, batched-vs-scalar local statistics, multi-query gene-level residual prediction, uncertainty, and exact preservation of observed coordinates.

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
- Removed prior-specific event fields and settings while preserving zero-only replacement, sparse output, and exact preservation of observed non-zero expression.
- Changed the default SuperCell-style membership graining level from `gamma = 20` to `gamma = 150` in both membership construction and the end-to-end workflow.

# DropoutKiller 0.2.0

- Rebuilt membership construction around the SuperCell graph/coarse-graining contract:
  biological hard strata -> embedding-space kNN graph -> walktrap hierarchy -> gamma-controlled memberships.
- Added membership-local ALRA-derived zero gating and negative-tail confidence scoring.
- Added event-only Gaussian neighbor borrowing and sparse selective replacement.
- Added Seurat integration, diagnostics, mathematical contract, unit tests, and GitHub Actions R CMD check.
- Enforced zero-only masks and exact preservation of observed non-zero expression.
