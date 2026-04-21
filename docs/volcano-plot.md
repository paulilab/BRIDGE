# Volcano Plot

Purpose:

- inspect differential expression for a selected contrast
- prioritize strongly changing and significant features

Inputs:

- one loaded dataset with processed output
- selected contrast
- p-value cutoff
- LFC cutoff
- optional feature highlight search

How to use:

1. Open `Volcano Plot`.
2. Choose contrast.
3. Set p-value and LFC cutoffs.
4. Optionally search/highlight features.
5. Click compute.

Outputs:

- interactive volcano plot
- significant-feature table for export

Interpretation tips:

- Upper-left / upper-right regions typically capture strong regulated features.
- Check that highlighted candidates behave as expected under chosen contrast.
- Use same thresholds when comparing across runs.

Suggested figure for this page:

- File: `docs/assets/figures/volcano_plot_example.png`
- Capture: computed volcano with at least one highlighted feature.
- Callouts: contrast selector, thresholds, highlighted label.
