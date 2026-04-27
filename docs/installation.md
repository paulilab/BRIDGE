# Installation

BRIDGE can be run with Docker or locally.

For the shortest path to a first plot, use this page together with [Quick Start](quick-start.md).

Before starting, you need a BRIDGE-compatible SQLite database. If you only have CSV/TSV files, check [Data Requirements](data-requirements.md) before creating the database with [Database Generation](database-generation.md).

## Fastest start (Docker)

Run:

```bash
docker run -d --rm \
  --name bridge \
  -p 3838:3838 \
  --mount type=bind,src=${YOUR_DATABASE},dst=/srv/data/database.db \
  ghcr.io/paulilab/bridge:latest
```

Then open:

```text
http://localhost:3838
```

## Local setup

Prerequisites:

- R
- Python
- A BRIDGE-compatible SQLite database

Clone the repository:

```bash
git clone https://github.com/paulilab/BRIDGE
cd BRIDGE
```

Restore the R environment:

```r
renv::restore()
```

Run the app:

```bash
Rscript app.R user_database.db
```

If you do not have a database yet, go to [Database Generation](database-generation.md).

## First-run success checklist

1. App opens at `http://localhost:3838`.
2. Species selector shows the expected species.
3. A data table appears after species selection.
4. A matching annotation table appears.
5. Dataset appears after clicking `Load Data`.
6. `PCA` or `Raw Heatmap` produces a plot.

After this works, load a second dataset and test `Raw Integration` if your goal is cross-omics comparison.
