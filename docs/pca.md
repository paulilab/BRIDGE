# PCA

Purpose:

- summarize global variance structure across samples
- inspect sample relationships and possible batch-like effects

Inputs:

- one loaded dataset with processed object

How to use:

1. Open `PCA`.
2. Click compute.
3. Inspect clustering/separation of samples.
4. Review loading table for feature contributions.

Outputs:

- PCA plot
- loading/contribution table

Interpretation tips:

- Clear group separation can indicate strong condition effects.
- Unexpected grouping can indicate technical variation or annotation issues.
- Pair PCA interpretation with raw/de views before filtering decisions.

Suggested figure for this page:

- File: `docs/assets/figures/pca_example_with_loadings.png`
- Capture: PCA plot and loadings table from the same dataset.
- Callouts: principal components and top contributing features.
