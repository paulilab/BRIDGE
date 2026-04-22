# Raw Heatmap

Purpose:

- visualize unprocessed expression patterns across selected datapoints
- detect broad trends and obvious quality issues

Inputs:

- one loaded dataset
- selected datapoint columns from `Loading Data`

How to use:

1. Open `Raw Heatmap` tab in a loaded dataset.
2. Click the compute button.
3. Inspect global intensity patterns across rows and columns.

What to look for:

- strong global shifts between stages/conditions
- outlier-like samples
- unusually flat or noisy signal regions

Notes:

- This module is for first-pass inspection.
- It does not apply significance filtering.

`<figure 1: Raw Heatmap view after compute, showing sample axis, feature axis, and overall signal pattern across selected datapoints.>`
