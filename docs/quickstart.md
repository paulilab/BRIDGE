# Quickstart

This is the fastest path to a working BRIDGE session.

Option A: run with Docker

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

Option B: run locally

```bash
git clone https://github.com/paulilab/BRIDGE
cd BRIDGE
```

```r
renv::restore()
```

```bash
Rscript app.R user_database.db
```

If you do not have a database yet:

```bash
touch user_database.db
python Python/db_adding.py
python Python/db_adding_annotation.py
```

First-session checklist:

1. Load one dataset in the Data Selection panel.
2. Open `Individual Exploration` and run one plot (for example Volcano Plot).
3. Load a second dataset and test `Raw Integration`.
