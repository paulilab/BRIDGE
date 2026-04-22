# Installation

BRIDGE can be run with Docker or locally.

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

## First-session checklist

1. Load one dataset in the Data Selection panel.
2. Open `Individual Exploration` and run one plot (for example Volcano Plot).
3. Load a second dataset and test `Raw Integration`.
