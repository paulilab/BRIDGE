# Raw Integration

Purpose:

- compare multiple loaded datasets at raw level
- inspect cross-omics trends before significance filtering

Requirements:

- at least two loaded datasets
- selected datapoint columns for each dataset
- equal number of selected datapoint columns across datasets

How to use:

1. Open `Raw Integration`.
2. Select datasets to integrate.
3. Choose datapoint columns for each dataset.
4. Check preview for column alignment.
5. Click integrate.

Outputs:

- integrated raw table (with source table label)
- integrated datapoints plot for selected genes

Practical tips:

- Use matching biological stages/conditions across datasets.
- Alignment is positional, so order matters.
- In phosphoproteomics, one gene can map to multiple peptide lines.

Suggested figure for this page:

- File: `docs/assets/figures/raw_integration_preview_and_plot.png`
- Capture: preview panel and integrated trend plot.
- Callouts: column matching, integrate action, combined output.
