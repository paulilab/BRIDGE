# Examples

This page gives practical user journeys for common BRIDGE tasks.

## Showcase dataset snapshot

If you use the showcase-style database, a common setup is:

| Table | Omics type | Rows | Datapoint columns |
| --- | --- | --- | --- |
| `zebrafish_proteomics_Showcase_1` | Proteomics | `273` | `18` |
| `zebrafish_phosphoproteomics_Showcase_1` | Phosphoproteomics | `533` | `18` |
| `zebrafish_rnaseq_Showcase_1` | RNA-seq | `1000` | `18` |

## Example 1: RNA-seq single-table exploration

Goal:

- inspect variance structure
- identify significant genes

Steps:

1. Load one RNA-seq table in `Loading Data`.
2. Open `PCA` and click compute.
3. Open `Volcano Plot`, choose contrast, set thresholds, compute.
4. Open `Gene Expression` and inspect selected genes over datapoints.

## Example 2: Proteomics table with enrichment

Goal:

- find significant signals and enriched pathways

Steps:

1. Load one proteomics table.
2. Run `DE Heatmap` and `Volcano Plot`.
3. Open `Enrichment Analysis`.
4. Pick contrast + database (GO, KEGG, or Reactome), set thresholds, compute.

## Example 3: Raw integration across datasets

Goal:

- compare trends across omics without significance filtering

Steps:

1. Load two or more datasets.
2. Open `Raw Integration`.
3. Select tables and matching datapoint columns.
4. Check integration preview.
5. Compute and inspect integrated raw table + datapoints plot.

## Example 4: Processed integration across matched contrasts

Goal:

- find intersecting significant IDs across datasets

Steps:

1. Load at least two datasets with computed processed outputs.
2. Open `Processed Integration`.
3. Select one contrast per dataset.
4. Set p-value/LFC thresholds and cluster number `k`.
5. Compute and inspect processed tables, heatmaps, and LFC scatterplots.

Tips:

- Use comparable contrasts across datasets.
- Start with moderate thresholds, then refine.
- Check preview dimensions before interpretation.
