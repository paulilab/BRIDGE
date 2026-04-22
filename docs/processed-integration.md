# Processed Integration

Purpose:

- integrate datasets based on significant differential signals
- focus on intersecting IDs across selected omics layers

Requirements:

- at least two loaded datasets
- one selected contrast per dataset
- p-value and LFC thresholds
- cluster number `k`

How to use:

1. Open `Processed Integration`.
2. Select datasets.
3. Select one comparison/contrast per dataset.
4. Set p-value and LFC thresholds.
5. Set `k`.
6. Click integrate.

Workflow summary:

1. Filter each dataset by selected contrast and thresholds.
2. Intersect significant IDs across datasets.
3. Build processed tables and matrices from intersected IDs.
4. Render integrated heatmaps and trend views.
5. Render LFC scatter comparisons where compatible.

Outputs:

- processed integrated tables
- integrated heatmaps
- cluster tables
- LFC scatterplot comparisons

Interpretation tips:

- Use biologically comparable contrasts across datasets.
- If intersection is too small, relax thresholds.
- Check preview dimensions before interpreting cluster patterns.

`<figure 8: Processed Integration workflow with settings panel, preview dimensions, integrated heatmaps, and LFC comparison output.>`
