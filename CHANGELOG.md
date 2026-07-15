# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.0] - 2026-07-14

First release: an ATAC-seq Snakemake workflow with **RPGC depth-normalized**
coverage (no spike-in). Adapted from the spike-in–normalized sibling workflow by
dropping spike-in alignment/normalization and using read-depth normalization
throughout.

### Added

- `workflow/Snakefile` — single standardized Snakemake Workflow Catalog entry
  point. A unified DAG builds the primary stage and the QC report in dependency
  order; run subsets with the `atacseq_all` and `qc_all` target rules.
- **Primary stage** (`atacseq_all`): FastQC/fastp → Bowtie2 alignment to a human
  index (`build_genome_index`) → filtering/dedup/blacklist → MACS2 peaks → RPGC
  depth-normalized bigWigs (`create_bigwig`) → reproducible fixed-width consensus
  peaks (majority vote / IDR by replicate count) with a featureCounts fragment
  matrix.
- **QC stage** (`qc_all`): deepTools coverage/fragment-size/fingerprint/
  correlation/PCA/GC/TSS, a numeric TSS-enrichment score, FRiP, IDR on relaxed
  peaks, library complexity (NRF/PBC1/PBC2), reads-in-annotation and peak
  summaries, a FastQC-only MultiQC report, and a self-contained interactive HTML
  QC report.
- **Differential analysis** (`ATACseq_Dx.ipynb`, R / Bioconductor): DESeq2
  (median-of-ratios) differential binding on the consensus counts, split into
  promoter vs distal peaks with a paired design, plus Gviz genome-browser tracks
  from the RPGC bigWigs.
- `config/` (config + sample sheet), `workflow/schemas/config.schema.yaml`
  (validated on every run and rendered as the catalog parameter table), `.test/`
  catalog test case, `workflow/envs/` per-rule conda environments, Docker image
  + `run_pipeline.sh`, CI, and `CITATION.cff`.

[Unreleased]: https://github.com/gynecoloji/snakemake_ATACseq/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/gynecoloji/snakemake_ATACseq/releases/tag/v1.0.0
