# DropoutKiller

`DropoutKiller` is an R package for **selective** scRNA-seq dropout recovery. It does not impute an entire expression matrix. It first identifies observed zeros with strong membership-local evidence of technical dropout, then modifies only those events.

The current biological-prior branch is deliberately restricted to **PPI and pathway information**. GRN / TF-target priors are not accepted by the public API.

## Complete workflow

```text
normalized scRNA-seq expression X (genes x cells)
                 |
                 +-----------------------------+
                 |                             |
                 v                             v
        low-dimensional embedding       broad biological strata
                 |                       (e.g. major cell type)
                 +-------------+---------------+
                               |
                               v
              SuperCell-style membership construction
       kNN graph -> walktrap hierarchy -> gamma-controlled cut
                               |
                               v
                   membership-local expression
                               |
                               v
                     local low-rank model
                               |
                               v
          ALRA-derived zero gate + negative-tail confidence
                               |
                    only observed X[g,c] = 0
                               |
                               v
                   high-confidence dropout mask
                               |
               +---------------+----------------+
               |                                |
               v                                v
      membership-local cell borrowing    PPI/pathway prior
        Gaussian latent-space weights      standardized locally
               |                                |
               +---------------+----------------+
                               |
                               v
                75% cell + 25% prior by default
                    (renormalized if one missing)
                               |
                               v
                    selective replacement only
                               |
                               v
                     recovered expression
```

## Design contract

1. **SuperCell-style membership**  
   Within trusted broad biological strata, cells are connected by a Euclidean kNN graph in a low-dimensional embedding. Walktrap hierarchical clustering is cut to approximately `n/gamma` memberships. Large strata can use anchor clustering plus nearest-centroid assignment.

2. **Membership-local ALRA gate**  
   Each membership is decomposed independently. A zero is eligible only if its low-rank reconstruction exceeds a gene-specific adaptive threshold estimated from the negative reconstruction tail.

3. **Confidence is evidence, not a posterior**  
   Negative reconstructed values define a working biological-zero null. The returned score is one-sided null-tail confidence and is not labelled as a calibrated posterior dropout probability.

4. **Cell-space recovery**  
   Candidate expression is borrowed only from cells in the same membership. Gaussian latent-space weights are used and, by default, only cells with observed positive expression for the target gene contribute.

5. **PPI/pathway prior only**  
   PPI and pathway adjacency matrices are aligned to the expression genes, self-edges are removed, and target rows are L1-normalized. Prediction is performed in membership-local standardized gene space to prevent high-abundance genes from dominating only because of scale.

6. **Selective hybrid recovery**  
   When both components exist: `Xhat = alpha * Xhat_cell + (1-alpha) * Xhat_prior`, with `alpha = 0.75` by default. If one component is unavailable, the remaining component is renormalized to weight 1.

7. **Observed values are immutable**  
   Any mask entry overlapping an observed non-zero value is rejected. The core workflow never overwrites observed non-zero expression.

## Installation

```r
remotes::install_github("1667857557/dropout-killer")
```

For the development PR branch:

```r
remotes::install_github(
  "1667857557/dropout-killer",
  ref = "feature/complete-selective-dropout-framework",
  upgrade = "never"
)
```

## 1. Matrix workflow

`x` must be a non-negative normalized expression matrix with genes in rows and cells in columns. `embedding` must contain cells in rows.

```r
library(DropoutKiller)

fit <- dropout_killer(
  x = x,
  embedding = pca,
  group = major_cell_type,
  split_by = condition,
  gamma = 20,
  k_knn = 5,
  rank = "auto",
  quantile_prob = 0.001,
  threshold = 0.95,
  alpha = 0.75,
  neighbor_k = 30,
  neighbor_positive_only = TRUE
)

fit
fit$expression
fit$events
validate_dropout_result(fit, x)
```

If condition-specific states must never borrow from one another, pass condition to `split_by`. If cross-condition borrowing within the same trusted biological state is intended, leave `split_by = NULL`.

## 2. Build PPI/pathway priors

### PPI

PPI edges are symmetric by default:

```r
ppi <- gene_prior_from_edges(
  ppi_edges,
  genes = rownames(x),
  source = "gene_a",
  target = "gene_b",
  weight = "score",
  prior_type = "ppi"
)
```

### Pathway

Pathway edges default to directed because some curated pathway resources encode source-to-target direction. Set `directed = FALSE` when the source does not provide interpretable directionality.

```r
pathway <- gene_prior_from_edges(
  pathway_edges,
  genes = rownames(x),
  source = "source_gene",
  target = "target_gene",
  weight = "weight",
  prior_type = "pathway"
)
```

The constructor accepts only `prior_type = "ppi"` or `"pathway"`.

## 3. Run with PPI only

```r
fit_ppi <- dropout_killer(
  x = x,
  embedding = pca,
  group = major_cell_type,
  split_by = condition,
  ppi = ppi,
  alpha = 0.75
)
```

## 4. Run with pathway only

```r
fit_pathway <- dropout_killer(
  x = x,
  embedding = pca,
  group = major_cell_type,
  split_by = condition,
  pathway = pathway,
  alpha = 0.75
)
```

## 5. Run with PPI + pathway

Both sources are aligned and row-normalized before fusion. By default they receive equal prior weight.

```r
fit_both <- dropout_killer(
  x = x,
  embedding = pca,
  group = major_cell_type,
  split_by = condition,
  ppi = ppi,
  pathway = pathway,
  prior_weights = c(ppi = 0.5, pathway = 0.5),
  alpha = 0.75
)
```

The `prior_weights` control only the composition of the 25% biological-prior branch. `alpha = 0.75` still controls cell-space versus biological-prior recovery.

## 6. Seurat workflow

```r
obj <- dropout_killer_seurat(
  object = obj,
  assay = "RNA",
  slot = "data",
  reduction = "pca",
  dims = 1:20,
  group_by = "major_cell_type",
  split_by = "condition",
  new_assay = "DropoutKiller",
  ppi = ppi,
  pathway = pathway,
  prior_weights = c(ppi = 0.5, pathway = 0.5)
)
```

Recovered values are continuous and are stored as **assay data**, not as raw counts.

## Main API

- `build_supercell_membership()`: graph-based membership construction.
- `local_alra_detect()`: membership-local low-rank zero detection.
- `select_dropout_mask()`: high-confidence zero-only sparse mask.
- `weighted_neighbor_prediction()`: membership-local cell prediction.
- `gene_prior_from_edges()`: construct a PPI or pathway adjacency prior.
- `prepare_gene_prior()`: align and normalize one PPI/pathway prior.
- `combine_gene_prior()`: fuse PPI and pathway priors.
- `gene_prior_prediction()`: event-level PPI/pathway prediction.
- `recover_dropout_expression()`: selective hybrid recovery for a supplied mask.
- `dropout_killer()`: end-to-end workflow.
- `dropout_killer_seurat()`: Seurat wrapper.
- `validate_dropout_result()`: validate core zero-preserving invariants.

See `inst/ALGORITHM.md` for the mathematical contract.

## Important trade-offs

- Hard `group`/`split_by` boundaries should use trusted broad labels. Over-fragmentation reduces membership size and local low-rank power.
- `threshold = 0.95` is a null-tail evidence threshold, not an FDR or calibrated posterior probability.
- Positive-only neighbor borrowing protects against unresolved technical zeros but estimates a conditional-positive local mean and can be upward-biased. Set `neighbor_positive_only = FALSE` for the literal all-neighbor Gaussian average.
- PPI/pathway information is **not used to decide whether a zero is a dropout**. It is evaluated only after the expression-derived dropout gate has passed.
- PPI/pathway priors can be context-mismatched. The prior branch therefore remains optional and subordinate to cell-space borrowing by default.
- Recovery is one-pass: recovered values are never fed back into membership construction or dropout detection.

## Methodological basis

The membership engine critically adapts the SuperCell graph/coarse-graining design: kNN graph, walktrap hierarchy, and gamma-controlled coarse graining. The zero gate and automatic-rank logic adapt ALRA's low-rank and negative-tail ideas. DropoutKiller does not copy either implementation verbatim; its explicit zero-only mask and PPI/pathway-only recovery contract define the intended behavior.
