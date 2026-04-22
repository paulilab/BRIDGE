# Enrichment Analysis

Purpose:

- interpret significant features through pathway/function context

Inputs:

- one loaded dataset with processed output
- selected contrast
- selected enrichment database
- p-value and LFC cutoffs for significant feature definition

Supported enrichment choices:

- GO
- KEGG
- Reactome

How to use:

1. Open `Enrichment Analysis`.
2. Choose contrast and enrichment database.
3. Set thresholds.
4. Click compute.

Outputs:

- enrichment dot plot

Important notes:

- Enrichment requires enough significant hits.
- If very few significant features pass thresholds, results may be empty.
- Phosphoproteomics enrichment may be unavailable depending on data/object context.

`<figure 5: Enrichment Analysis controls and resulting dot plot, showing contrast choice, enrichment database, and significant-term output.>`
