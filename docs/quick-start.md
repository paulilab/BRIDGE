# Quick Start

This page gets you from opening BRIDGE to generating one successful plot.

Goal:

- start the app
- load one dataset
- confirm the dataset is available
- generate a first plot

If you want the conceptual background before running anything, read [Introduction](intro.md) and [Methodology](methodology.md) first.

## Choose your route

| Starting point | Go here first |
| --- | --- |
| You already have a BRIDGE-compatible SQLite database | Continue with this page |
| You have CSV/TSV files but no database | Read [Data Requirements](data-requirements.md), then [Database Generation](database-generation.md) |
| You were given a showcase database | Use that database for the steps below |
| You want the exact processing/integration logic | Read [Methodology](methodology.md) |

## Step 1: start BRIDGE

Docker:

```bash
docker run -d --rm \
  --name bridge \
  -p 3838:3838 \
  --mount type=bind,src=${YOUR_DATABASE},dst=/srv/data/database.db \
  ghcr.io/paulilab/bridge:latest
```

Local:

```bash
Rscript app.R user_database.db
```

Then open:

```text
http://localhost:3838
```

Expected result:

- the BRIDGE interface opens in your browser
- the Data Selection panel is visible

If the app does not open, return to [Installation](installation.md).

## Step 2: load one dataset

In the Data Selection panel:

1. Select the species.
2. Select one data table.
3. Select datapoint columns.
4. Select the matching annotation table.
5. Click `Load Data`.

Expected result:

- the dataset appears in the loaded-datasets area
- the individual exploration modules become available
- no table is deleted from the database when it is removed from the session

If no table or annotation appears, check [Loading Data](loading-data.md) and [FAQ](faq.md).

## Step 3: generate your first plot

Recommended first plot:

1. Open `PCA`.
2. Click compute.
3. Confirm that a PCA plot appears.

Alternative first plot:

1. Open `Raw Heatmap`.
2. Click compute.
3. Confirm that a heatmap appears.

Why start here:

- PCA and Raw Heatmap give a global view of the loaded dataset.
- They are good first checks before interpreting differential results.
- They can reveal obvious sample grouping, outliers, or data-loading issues.

## Step 4: decide what to do next

| Goal | Next module |
| --- | --- |
| Check broad sample structure | [PCA](pca.md) |
| Inspect raw expression patterns | [Raw Heatmap](raw-heatmap.md) |
| Find differential signals | [Volcano Plot](volcano-plot.md) or [DE Heatmap](de-heatmap.md) |
| Inspect candidate genes/features | [Gene Expression](gene-expression.md) |
| Interpret significant features biologically | [Enrichment Analysis](enrichment-analysis.md) |
| Compare raw trends across datasets | [Raw Integration](raw-integration.md) |
| Compare shared significant IDs across datasets | [Processed Integration](processed-integration.md) |

## First-run success checklist

- App opens at `http://localhost:3838`.
- The species selector shows the expected species.
- A data table appears after species selection.
- A matching annotation table appears.
- A dataset appears after clicking `Load Data`.
- `PCA` or `Raw Heatmap` produces a plot.

## If the first run fails

Start with these checks:

- Docker users: confirm `${YOUR_DATABASE}` points to the intended SQLite file.
- Confirm the data table name starts with the selected species.
- Confirm the table is registered in `table_metadata`.
- Confirm the annotation table is registered in `annotation_metadata`.
- Confirm identifier columns are present and non-empty.
- Confirm datapoint columns follow the replicate suffix pattern, for example `X6.hpf_1`.

For formatting rules, see [Data Requirements](data-requirements.md).
