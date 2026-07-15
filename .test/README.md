# Catalog test case (`.test/`)

This directory exists so the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/)
can render the workflow's **tube map** (rule graph) with
[snakevision](https://github.com/snakemake/snakevision). The catalog runs:

```bash
snakemake -s workflow/Snakefile -c 1 -d .test --forceall --rulegraph
```

`--rulegraph` only resolves the **rule dependency graph** — no rule is executed —
so the reference genomes and FASTQs here are **empty 0-byte placeholders**. They
exist only so the DAG resolves; nothing reads their contents. The sample sheet
uses two 2-replicate groups so the IDR-dependent rules also appear in the map.

Contents:

- `config/config.yaml` — copy of the top-level config (paths resolve under `.test/`)
- `config/samples.csv` — 4 samples, 2 groups × 2 replicates
- `data/*.fastq.gz` — empty placeholder paired-end reads
- `ref/*` — empty placeholder reference genomes / annotations / BEDs

**This is not an end-to-end integration test.** Turning it into one (so the
workflow also earns the catalog's "tests" ranking) would require real miniature
reference genomes + reads that actually run through the tools, plus `picard.jar`
and the per-rule conda environments — a larger, separate effort.
