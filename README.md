# DropoutKiller

DropoutKiller detects recoverable RNA zeros and selectively imputes those calls.
Version 0.8 defaults to `Supercell_hierarchy_Lean_membership`. Paired RNA+ATAC
Seurat objects automatically use `Supercell_hierarchy_Lean_membership_WNN`.
Observed nonzero expression remains unchanged on the library-normalized log scale.

## Usage

```r
# broad_type contains major lineages, NOT fine cell subclusters.
# Raw RNA counts are used even when ATAC is the object's active assay.
out <- dropout_killer_seurat(object, group_by = "broad_type",
                             return_result = TRUE)
out$result$settings$detection_method
model <- out$result$detection$model

# Reuse a separately validated model of the same modality.
out2 <- dropout_killer_seurat(object, group_by = "broad_type",
                              lean_model = model, return_result = TRUE)
```

For RNA matrices, supply raw counts, a cell-by-dimension RNA embedding and broad
groups to `dropout_killer(counts, embedding, group = broad_group)`. The supplied
embedding controls recovery; the default RNA detector always rebuilds its geometry
from the normalized RNA SVD, matching the thinning calibration pipeline. This also
applies when reusing a calibrated model. Explicit custom `lean_geometry` requires
an externally fitted model calibrated with that same geometry builder. In the Seurat
RNA route, the recovery embedding is also rebuilt from RNA SVD. Fine cell labels
are never selected automatically. `group_by`/`group` must be explicitly chosen;
with a supplied calibration model, omitting groups treats all cells as one group.
The package cannot infer which annotation represents a broad biological lineage.

## Detection and hierarchy

1. Library-size normalize RNA to 10,000 and apply log1p. Retain the Lean benchmark's
   ALRA singular-value-spacing rank rule and gene-wise standardized margin `z`.
2. Build the neighborhood graph within each supplied broad class, retain the
   complete Walktrap merge history, and cut at `floor(n / gamma)` memberships
   (`gamma=150`). The cut is bounded below by one and the connected-component
   count. No cross-class merge or invented connection is introduced.
3. Select up to `neighbor_k=30` support cells from the query's own membership.
   Only if support is insufficient, borrow from the closest memberships in the
   retained hierarchy. RNA ties use centroid distance, then stable IDs; WNN ties
   use affinity where available, then stable IDs. Self-neighbors are excluded.
4. Compute weighted positive-neighbor support and the benchmark's smoothed logit
   `(support + 0.01) / 1.02`. No available neighbor contributes neutral support 0.5.
5. Compute the membership term
   `H = n_m/(n_m+50) * (logit((positive_m+0.5)/(n_m+1)) -
   logit((positive_group+0.5)/(n_group+1)))`.
6. Fit or reuse a class-balanced logistic model on `z`, `local_hierarchy`, and `H`.
   Call only zeros whose linear score is strictly above the calibrated threshold.
   The sigmoid score is not a calibrated biological-dropout posterior.

The P1 stabilized-state recovery engine remains the default recovery stage.
`threshold` continues to govern historical confidence detectors; Lean uses the
stored model threshold and its `fpr` calibration target. Returned detection
events and optional sparse scores cover called zeros, not every uncalled zero.

## Paired RNA and ATAC

The Seurat wrapper detects an unambiguous ChromatinAssay or assay named ATAC/peaks.
Specify `atac_assay` if multiple ATAC assays exist. RNA+ADT is not automatically
treated as RNA+ATAC. `modality="rna"` explicitly opts out of WNN. RNA-only and
explicit historical routes do not inspect or validate unused ATAC assays.

Every WNN build recomputes RNA NormalizeData/variable features/ScaleData/PCA and
ATAC TF-IDF/top features/LSI. It uses PCA 1:40 and LSI 2:40, reduced to available
dimensions for small inputs. LSI component one is always excluded. Seurat's
FindMultiModalNeighbors is rerun within each broad class; old graphs and reductions
are not reused. Defaults are `wnn_npcs=40`, `wnn_k=30`.

The graph follows the actual SuperCell 2.0 `ComputeMultimodalKnn` kernel:
`A[i,j] = 1 - 2 * weighted.nn@nn.dist[i,j]^2`, followed by `A + t(A)`, removing
self-loops and clamping negative roundoff weights to 1e-16 before Walktrap.
It does not substitute the shared-neighbor `wsnn` graph or concatenated PCA/LSI.
Hierarchy support uses positive WNN affinities and row-normalizes selected edges;
it does not fill absent WNN edges with RNA Euclidean neighbors. Strata with fewer
than four cells are explicitly warned about and retain singleton memberships and
neutral neighbor support. Seurat/Signac failures propagate; no silent modality
fallback occurs. RNA and ATAC must have identical paired cell sets.

## Calibration

If `lean_model` is absent, default calibration uses q=0.5 UMI retention and three
independent seeds (10001:10003), rebuilding geometry for every thinned RNA matrix.
It trains on all artificially lost positives plus all eligible proxy negatives;
there is no event subsampling. Proxy negatives are original zeros whose gene has
prevalence <=0.005 in their class and >=0.20 in another broad class. Calibration
requires raw integer counts and at least two broad classes. It stops if suitable
negatives do not exist. Any relaxation must be explicit, for example
`lean_control=list(negative_max=0.01)`; do not reinterpret it as biological truth.

`lean_control` accepts `calibrate_lean_detector` options such as `q`, `seeds`,
`fpr`, `negative_max`, and `negative_other_min`. Default `fpr=0.01` is an empirical
training proxy-FPR, not a guarantee on unseen data or a biological false-discovery
rate. Use separate masks/data for evaluation. Calibration can be expensive for
large objects, especially WNN; save and reuse an independently validated model.
For normalized inputs, supply `lean_model` explicitly. Fitted coefficients are
not silently imported from a PBMC/TNBC benchmark fold or shared between modalities.

The lower-level `fit_lean_detector(features, truth, method=...)` supports external
training with the same three feature definitions. `supercell_lean_detect` runs
detection alone. `build_wnn_supercell` returns inspectable WNN geometry/provenance.

## Compatibility and evidence

Historical detectors remain explicit options: `alra_global_by_group`,
`eb_zero_null`, and `alra_quantile`. The new detector builds its own retained
hierarchy: a bare membership vector or `split_by` is rejected, rather than losing
the broad-class-only contract. Existing legacy tests select their old detector
explicitly; separate tests exercise the new default and automatic WNN route.

The earlier local WNN benchmark used a wsnn graph and broad-group edges without
strict shortage-based hierarchy borrowing. Its scores cannot validate this
corrected kernel/hierarchy implementation. RNA graph cuts also now follow the
SuperCell floor rule rather than the older round rule. Full dataset performance
must be remeasured before claiming equivalence to those archived tables.

Actual upstream source inspected for this implementation:
- [ComputeMultimodalKnn](https://github.com/GfellerLab/SuperCell/blob/89b34c078dba0b91289a5a23cf75a6db3a46d232/R/SuperCell_for_Seurat.R)
- [SCimplify_from_Seurat / Walktrap and floor cut](https://github.com/GfellerLab/SuperCell/blob/89b34c078dba0b91289a5a23cf75a6db3a46d232/R/SCimplify_for_Seurat_v5.R)

This is an implementation of the graph and hierarchy conventions, not a dependency
on the SuperCell package. Signac and Seurat are optional dependencies needed only
for the multiome wrapper. Assay5 split count layers must first be joined explicitly.
