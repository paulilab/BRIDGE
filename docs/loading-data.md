# Loading Data

The Data Selection panel controls what becomes available in BRIDGE.

Main controls:

1. Select species.
2. Select data table.
3. Select datapoint columns.
4. Select annotation table.
5. Click `Load Data`.
6. Optionally remove loaded tables.

How BRIDGE uses your selections:

- It reads table and annotation metadata from SQLite.
- It queries selected columns from the chosen data table.
- It computes or loads cached processed output.
- It enables module tabs for the loaded dataset.

Loaded datasets:

- are listed in the sidebar
- can coexist in one session
- can be removed from the session without deleting from the database

Common issues:

- No tables appear:
  check table naming, species prefix, and `table_metadata`.
- No annotation appears:
  check `annotation_metadata`.
- Load succeeds but plots are empty:
  validate identifiers and selected datapoints.

First successful load checklist:

- The selected species matches the table-name prefix.
- A data table is visible in the table selector.
- Datapoint columns are selected in the intended biological order.
- A matching annotation table is selected.
- The loaded dataset appears in the session after clicking `Load Data`.
- `PCA` or `Raw Heatmap` can compute a first plot.

Practical tip:

- Start with one table and a small datapoint subset to validate setup quickly.
