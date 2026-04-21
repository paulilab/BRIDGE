# Introduction

BRIDGE is built for users who want to analyze and integrate omics datasets through an interface instead of writing analysis scripts for every step.

Core idea:

- Data live in a local SQLite database.
- You load selected tables and columns through the sidebar.
- BRIDGE creates per-dataset analysis workspaces.
- You can compare datasets through raw and processed integration.

Two main workflows:

- Individual Exploration:
  inspect one loaded dataset at a time with module tabs.
- Integration:
  combine multiple loaded datasets to compare trends and significant signals across omics layers.

Typical session:

1. Generate or open a BRIDGE-compatible database.
2. Load one or more datasets.
3. Explore each dataset individually.
4. Run raw or processed integration.
5. Export tables/plots for reporting.

Scope note:

- The backend is species-aware through naming conventions.
- The current UI species selector is configured around zebrafish workflows.
