# Installation

BRIDGE can be run with Docker or locally.

Docker run:

```bash
docker run -d --rm \
  --name bridge \
  -p 3838:3838 \
  --mount type=bind,src=${YOUR_DATABASE},dst=/srv/data/database.db \
  ghcr.io/paulilab/bridge:latest
```

Local run prerequisites:

- R
- Python
- A BRIDGE-compatible SQLite database

Local setup:

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

Optional documentation setup:

```bash
python -m pip install -r requirements-docs.txt
python -m mkdocs serve
```

If `mkdocs serve` uses the wrong binary, run:

```bash
python -m mkdocs serve
```
