# FAQ

## Which omics types are supported?

Proteomics, phosphoproteomics, and RNA-seq.

## Why does my table not appear after selecting species?

Most common reasons:

- table name does not follow expected species prefix
- `table_metadata` entry is missing
- table naming pattern is inconsistent with the selected species

## Which identifier columns are required?

For data tables:

- `Gene_Name`
- `Gene_ID`
- `Protein_ID`

For phosphoproteomics, add:

- `pepG`

## Why is annotation required?

Annotation is used for gene-level context and transformations. BRIDGE expects an annotation table that matches species and genome version.

## What is the expected datapoint naming format?

Datapoint columns should end with replicate index using underscore, for example:

```text
X6.hpf_1
```

## Why do no genes appear after filtering?

Likely causes:

- thresholds too strict
- selected contrast has weak signal
- poor ID overlap across integrated datasets

Try more permissive thresholds first, then tighten.

## What is the difference between raw and processed integration?

- Raw integration compares loaded values directly across selected columns.
- Processed integration filters by significance first, then intersects IDs.

## Can I load multiple datasets in one session?

Yes. This is required for integration workflows.

## Can I use preprocessed objects?

Yes. During database generation you can add precomputed processed objects (for example `SummarizedExperiment`) so BRIDGE can reuse them.
