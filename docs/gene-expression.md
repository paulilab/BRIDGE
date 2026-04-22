# Gene Expression

Purpose:

- inspect expression trajectories for selected genes across datapoints
- compare one or multiple genes in the same view

Inputs:

- one loaded dataset
- one or more selected genes
- selected transformation scale

Available transformations:

- Proteomics/phosphoproteomics:
  continuous, log-scale, total intensity, median normalization
- RNA-seq:
  continuous, log-scale, TPM, FPKM, TMM, CPM

How to use:

1. Open `Gene Expression`.
2. Search and select one or more genes.
3. Choose transformation/scale.
4. Inspect line/trend behavior across stages.

Outputs:

- gene-level trend plot
- corresponding table output for selected view

Practical tips:

- For phosphoproteomics, one gene may produce multiple peptide trends.
- Validate stage ordering and replicate naming in your input columns.
- Compare transformed and continuous views before conclusions.

`<figure 4: Gene Expression module with multiple selected genes, visible transformation selector, and stage-wise trend lines.>`
