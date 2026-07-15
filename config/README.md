# Configuration

This workflow is configured through two files in this directory:

- `config.yaml` — all workflow parameters (see below)
- `samples.csv` — the sample sheet

plus reference data you download into `ref/` (not tracked in git; see
[Reference data](#reference-data)).

## Sample sheet (`config/samples.csv`)

CSV with one row per sample and these columns:

| column      | description                                                                 |
|-------------|-----------------------------------------------------------------------------|
| `sample_id` | Sample name. Raw reads must be `data/<sample_id>_R1_001.fastq.gz` / `_R2_001.fastq.gz`. |
| `type`      | Free-text condition label (e.g. `Control`, `NICD3`).                         |
| `group`     | Replicate group. Reproducibility handling is chosen from group size (below). |

Example:

```csv
sample_id,type,group
GSF4007-Control_1_S11,Control,group1
GSF4007-Control_2_S13,Control,group1
GSF4007-Control_3_S15,Control,group1
GSF4007-NICD3-V5_1_S12,NICD3,group2
GSF4007-NICD3-V5_2_S14,NICD3,group2
GSF4007-NICD3-V5_3_S16,NICD3,group2
```

**Per-group reproducibility** is derived automatically from the number of
replicates in each `group`:

- **≥ 3 replicates** → majority vote (a peak is kept if it recurs in ≥
  `consensus_min_replicates` replicates).
- **exactly 2 replicates** → IDR (`idr_threshold`).
- **1 replicate** → the sample's own peaks are used as-is.

## Parameters (`config/config.yaml`)

Every parameter — with its type, default, and description — is defined once in the
config schema, [`workflow/schemas/config.schema.yaml`](../workflow/schemas/config.schema.yaml).
That schema is the single source of truth: the workflow validates `config.yaml`
against it on every run (and fills in defaults for anything you omit), and the
[Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=gynecoloji/snakemake_ATACseq)
renders it as a parameter table on the workflow page.

To configure a run, edit `config.yaml` directly — it ships with working defaults
and an inline comment on every parameter. At minimum, point the reference-file
paths (`human_fasta`, `blacklist`, `gtf`, `promoter_bed`, `enhancer_bed`) at the
files you provide (see [Reference data](#reference-data)).

## Reference data

Genomes, indexes and large annotations are **not** shipped in the repo (they are
`.gitignore`d). Download / place them under `ref/` before running, matching the
paths in `config.yaml`:

- `ref/hg38.fa` — chr-prefixed UCSC human genome
- `ref/hg38_blacklist_regions.bed` — ENCODE hg38 blacklist (shipped)
- `ref/gencode.v36.annotation.gtf` — GENCODE annotation (for TSS QC)
- `ref/hg38.2bit` — for `computeGCBias`
- `ref/picard.jar` — Picard (used by MarkDuplicates)

The Bowtie2 index (`ref/BOWTIE2/`) is built automatically by the
`build_genome_index` rule from `human_fasta`.

See the top-level `README.md` for full setup and run instructions.
