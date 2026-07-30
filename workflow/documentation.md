# Technical documentation

Step-by-step documentation of the `snakemake_ATACseq` workflow. For
installation, container usage, and the full narrative, see the top-level
[`README.md`](../README.md); for every configuration parameter, see
[`config/README.md`](../config/README.md) and the schema
[`workflow/schemas/config.schema.yaml`](schemas/config.schema.yaml).

The rule graph is shown in [`images/rulegraph.svg`](../images/rulegraph.svg) and
rendered as a "tube map" on the workflow's Snakemake Workflow Catalog page.

## Overview

A single `snakemake -s workflow/Snakefile --use-conda` run builds a **unified
DAG** covering two stages in dependency order:

1. **Primary** (`atacseq_all` target) — alignment → filtering → peak calling →
   RPGC depth-normalized coverage → reproducible consensus peaks + fragment counts.
2. **QC** (`qc_all` target) — deepTools QC, FRiP, IDR, library complexity,
   TSS enrichment, and an interactive HTML QC report.

Two **optional** stages are **not** in the default DAG — run them on demand after
the primary stage: **differential openness** (`diffopen_all`,
`workflow/rules/diffopen.smk`) and **TF footprinting** (`footprint_all`,
`workflow/rules/footprint.smk`). See the per-stage step sections below.

## Inputs

| Input | Location | Notes |
|---|---|---|
| Paired-end reads | `data/<sample_id>_R1_001.fastq.gz`, `_R2_001.fastq.gz` | one pair per sample |
| Sample sheet | `config/samples.csv` | columns `sample_id, type, group` |
| Human genome FASTA | `ref/hg38.fa` | chr-prefixed UCSC |
| Blacklist BED | `ref/hg38_blacklist_regions.bed` | ENCODE, chr-prefixed |
| GTF / 2bit / promoter+enhancer BEDs | `ref/…` | QC references |
| Picard | `ref/picard.jar` | duplicate marking |

Configuration is read from `config/config.yaml` and validated against the schema
at parse time (missing/invalid parameters fail fast).

## Steps (primary stage)

1. **`fastqc`** — raw-read quality.
2. **`fastp`** — adapter trimming + quality filtering (auto-detects adapters).
3. **`build_genome_index`** — optionally subset the human genome to the requested
   chromosomes (`align_chroms`) and build the Bowtie2 index.
4. **`bowtie2_align`** — alignment to the human index.
5. **`samtools_sort_filter_index`** — keep uniquely-mapped, properly-paired
   reads; record mitochondrial-% QC; restrict to the analysis chromosomes.
6. **`remove_duplicates`** — Picard MarkDuplicates.
7. **`filter_blacklist`** — fragment-level ENCODE blacklist removal.
8. **`call_peaks`** — MACS2 (BAMPE, `-q 0.05`).
9. **`create_bigwig`** — RPGC depth-normalized coverage bigWig.
10. **Consensus peaks** (`relaxed_peaks`, `reproducible_idr`, `consensus_peaks`,
    `count_fragments_consensus`) — per-group reproducibility (majority vote for
    ≥3 replicates, IDR for exactly 2), a fixed-width consensus set, and a
    featureCounts fragment matrix.

## Steps (QC stage)

deepTools coverage/fragment-size/fingerprint/correlation/PCA/GC/TSS, a numeric
TSS-enrichment score, FRiP, IDR on relaxed peaks, library complexity
(NRF/PBC1/PBC2), reads-in-annotation and peak summaries, a FastQC-only MultiQC
report, and a self-contained interactive HTML QC report
(`results/qc/atacseq_qc_report.html`).

## Steps (differential-openness stage, optional)

Opt-in DESeq2 differential chromatin openness (`diffopen_all` target) on the
consensus count matrix; contrast from the sample sheet's `type` column (reference
= `diffopen_ref_label`). Runs under spike-in-free normalization modes
(`diffopen_modes`, default `none` + `ctcf`), each in its own
`results/diffopen/<mode>/` directory. Because there is no spike-in, every mode
treats its anchor set as invariant and cannot detect a uniform global shift —
read the results as relative and compare modes.

1. **`diffopen`** (per mode) — DESeq2 size factors (`none` = median-of-ratios over
   all peaks; `ctcf` = restricted to constitutive CTCF anchors), promoter/enhancer/
   all split, MA plot, nominal p05/p01 subsets, `run_summary.txt`.
2. **`diffopen_annotate`** — nearest-TSS gene assignment (ChIPseeker) per mode.
3. **`diffopen_enrich`** — offline GO enrichment (clusterProfiler) per mode.
4. **`diffopen_bigwig`** + **`diffopen_tracks`** — per-mode size-factor-scaled
   bigWigs and Gviz browser tracks for the top regions.
5. **`diffopen_report`** — `results/diffopen/diffopen_report.html` comparing the
   normalizations side by side.

Run with: `snakemake -s workflow/Snakefile --use-conda --cores N diffopen_all`.

## Steps (footprinting stage, optional)

Opt-in TOBIAS differential TF footprinting (`footprint_all` target); requires a
JASPAR motif file at `config['jaspar_motifs']`. Conditions are the sample
sheet's `type` values.

1. **`footprint_merge_condition`** — `samtools merge` the blacklist-filtered BAMs
   of each condition into one pooled BAM.
2. **`tobias_atacorrect`** — model + remove Tn5 insertion bias over the consensus
   peaks → bias-corrected signal bigWig per condition.
3. **`tobias_scorebigwig`** — footprint-score track per condition.
4. **`tobias_bindetect`** — scan JASPAR motifs across conditions → differential
   TF-binding table, volcano figures, and per-TF footprints/BEDs
   (`results/footprint/bindetect/`).

Run with: `snakemake -s workflow/Snakefile --use-conda --cores N footprint_all`.

## Outputs

All outputs are written under `results/` (peaks, bigWigs, consensus matrix, QC
tables and reports); per-rule logs under `logs/`. See the README's "Output
Structure" section for the full tree.

## Running the tests

```bash
python -m pytest tests/ -q                               # unit tests
snakemake -s workflow/Snakefile -c 1 -d .test --forceall --rulegraph   # DAG/tube map
```
