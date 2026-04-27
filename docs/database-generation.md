# Database Generation

BRIDGE uses a SQLite database as backend. This database stores:

- raw tables
- annotation tables
- metadata tables
- optional cached processed objects

![BRIDGE SQLite backend structure](assets/sqlite_backend_structure.svg)

*Database figure. The database-generation scripts upload raw data and annotation tables into SQLite, then register metadata tables that BRIDGE uses to discover available datasets, annotation versions, and optional cached processed objects.*

<details>
<summary>Help: how to read this database diagram</summary>
<p>The raw and annotation tables store the biological data. The metadata tables tell BRIDGE which tables exist, what species or datatype they belong to, and how they should be loaded. If a table exists in SQLite but is missing from the metadata, it may not appear in the app.</p>
</details>

Create a new database file:

```bash
touch user_database.db
```

Step 1: add data tables

```bash
python Python/db_adding.py
```

What this script does:

1. Loads a CSV/TSV file.
2. Lets you select identifier columns.
3. Lets you select datapoint columns.
4. Uploads the table to SQLite.
5. Registers metadata in `table_metadata`.
6. Optionally stores a processed R object in cache.

Step 2: add annotation tables

```bash
python Python/db_adding_annotation.py
```

What this script does:

1. Loads annotation CSV/TSV.
2. Uploads annotation table.
3. Registers it in `annotation_metadata`.

Re-running scripts:

- You can run both scripts as many times as needed.
- Uploading with the same table name replaces that table.
- Metadata entries are updated with `INSERT OR REPLACE`.

Optional CLI usage:

`db_adding.py` supports non-interactive flags such as:

- `--csv`
- `--db`
- `--table`
- `--id-cols`
- `--tp-cols`
- `--processed`
- `--rds`

`db_adding_annotation.py` supports:

- `--csv`
- `--db`
- `--table`
