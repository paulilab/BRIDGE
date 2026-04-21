# DE Heatmap

Purpose:

- visualize differentially expressed features after processing
- explore cluster structure and significant patterns

Inputs:

- one loaded dataset
- selected contrast context from processed results
- p-value and LFC thresholds
- optional clustering settings and `k`

How to use:

1. Open `DE Heatmap`.
2. Set p-value and LFC cutoffs.
3. Enable/adjust clustering if needed.
4. Click compute.

Outputs:

- heatmap of significant features
- cluster-oriented view depending on settings
- table output for inspected values

Practical tips:

- Start with moderate thresholds to avoid empty output.
- If clusters are unstable, try fewer/more conservative features first.
- Use alongside Volcano Plot for consistency checks.

Suggested figure for this page:

- File: `docs/assets/figures/de_heatmap_thresholds.png`
- Capture: DE heatmap settings + rendered output.
- Callouts: thresholds, clustering toggle, resulting cluster blocks.
