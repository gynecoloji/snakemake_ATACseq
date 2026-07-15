[![CI](https://github.com/gynecoloji/snakemake_ATACseq/actions/workflows/ci.yml/badge.svg)](https://github.com/gynecoloji/snakemake_ATACseq/actions/workflows/ci.yml)

## Citation

If you use this workflow in your research, please cite it. Use the **"Cite this
repository"** button on the GitHub repository page (generated from
[`CITATION.cff`](CITATION.cff)). A DOI will be minted for the archived release on
the first tagged release.

**Please also cite the individual tools used:**
- **Snakemake**: Köster, J. and Rahmann, S. (2012). Snakemake—a scalable bioinformatics workflow engine. Bioinformatics, 28(19), 2520-2522.
- **MACS2**: Zhang, Y. et al. (2008). Model-based analysis of ChIP-Seq (MACS). Genome Biology, 9, R137.
- **deepTools**: Ramírez, F. et al. (2016). deepTools2: a next generation web server for deep-sequencing data analysis. Nucleic Acids Research, 44(W1), W160-W165.
- **Bowtie2**: Langmead, B. and Salzberg, S.L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9, 357-359.
- **SAMtools**: Li, H. et al. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078-2079.
- **IDR**: Li, Q. et al. (2011). Measuring reproducibility of high-throughput experiments. Annals of Applied Statistics, 5(3), 1752-1779.
- **featureCounts (Subread)**: Liao, Y. et al. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7), 923-930.
- **Consensus peaks (fixed-width / score-per-million)**: Corces, M.R. et al. (2018). The chromatin accessibility landscape of primary human cancers. Science, 362(6413), eaav1898.
- **BEDTools**: Quinlan, A.R. and Hall, I.M. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics, 26(6), 841-842.
- **FastQC**: Andrews, S. (2010). FastQC: a quality control tool for high throughput sequence data.
- **fastp**: Chen, S. et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884-i890.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

**Author**: gynecoloji  
**Project Repository**: [https://github.com/gynecoloji/snakemake_ATACseq](https://github.com/gynecoloji/snakemake_ATACseq)

For questions, issues, or feature requests, please:
1. Check the existing [Issues](https://github.com/gynecoloji/snakemake_ATACseq/issues) on GitHub
2. Submit a new issue with detailed information about your problem
3. Include relevant log files and system information for troubleshooting

## Acknowledgments

This pipeline was developed based on best practices from the ENCODE consortium and incorporates methodologies from multiple published ATAC-seq analysis workflows. Special thanks to the developers of all the integrated tools that make this comprehensive analysis possible.

---

**Note**: This pipeline is optimized for human genome analysis (hg38) but can be adapted for other organisms by updating reference files and parameters accordingly.# ATAC-seq Analysis Pipeline

A comprehensive Snakemake workflow for processing and analyzing ATAC-seq data from raw reads to peak calling with extensive quality control metrics and differential binding analysis.

## Overview

This pipeline integrates three complementary components for complete ATAC-seq analysis:

1. **Primary ATAC-seq stage** (`atacseq_all` target) - Raw FASTQ → Bowtie2 alignment to the human genome → filtering → MACS2 peak calling → **depth-normalized (RPGC) coverage tracks** → a reproducible, fixed-width **consensus peak set** with a fragment-count matrix
2. **ATAC-seq QC stage** (`qc_all` target) - deepTools QC, FRiP, IDR, library complexity, TSS enrichment score, a self-contained **interactive HTML QC report** (all QC except FastQC), plus a FastQC-only MultiQC

Both stages live in a single standard-layout `workflow/Snakefile`: one `snakemake --use-conda` run builds the primary stage **and** the QC report in dependency order (unified DAG). Run a subset with the `atacseq_all` or `qc_all` targets. The layout follows the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/) conventions, so the workflow can be deployed into another project with `snakedeploy deploy-workflow` (see [Deploying with snakedeploy](#deploying-with-snakedeploy)).
3. **Differential Analysis Notebook** (`ATACseq_Dx.ipynb`, R / Bioconductor) - **DESeq2** differential binding (NICD3 vs Ctrl) on the consensus count matrix — split into **promoter vs distal** peaks with a **paired** design under **default (median-of-ratios) normalization** — plus **Gviz genome-browser tracks** at Notch target loci as a positive control; runs on the `atacseq-diffbind` container

## Workflow Diagram

The workflow rule graph, rendered as a "tube map" with
[snakevision](https://github.com/OpenOmics/snakevision):

![ATAC-seq workflow tube map](images/rulegraph.svg)

The same tube map is rendered automatically on the
[Snakemake Workflow Catalog page](https://snakemake.github.io/snakemake-workflow-catalog/?usage=gynecoloji/snakemake_ATACseq)
from the executable test case in [`.test/`](.test). Regenerate it with:

```bash
snakemake -s workflow/Snakefile -c 1 -d .test --forceall --rulegraph > rulegraph.dot
snakevision -s all atacseq_all qc_all -o images/rulegraph.svg rulegraph.dot
```

## Features

- **Complete end-to-end processing** of paired-end ATAC-seq data
- **Depth-normalized coverage** (RPGC / 1×, "reads per genomic content") bigWig tracks for cross-sample comparison
- **Bowtie2 alignment** to the human genome
- **Chromosome control**: align to a configurable set (default chr1–22, chrX, chrM), record mitochondrial-% QC, then keep only the analysis chromosomes
- **Reproducible consensus peaks** (Corces-2018 fixed-width / score-per-million; majority-vote or IDR reproducibility) + a **featureCounts** fragment matrix for differential analysis
- **Blacklist region filtering** for removal of technical artifacts
- **Extensive QC metrics** (fragment size, TSS heatmap + numeric enrichment score, fingerprint, correlation/PCA, GC bias, FRiP, NRF/PBC1/PBC2, reads in promoters/enhancers)
- **Interactive HTML QC report** (self-contained, theme-aware; all QC except FastQC) + a FastQC-only MultiQC
- **Signal track generation** (bigWig, bedGraph) for visualization
- **Conda environment management** (one env per rule) plus a ready-to-run **Docker / Apptainer** image

## Pipeline Components

### 1. Primary Processing Pipeline (`atacseq_all` target)

**Processing Steps:**
```
Raw FASTQ → FastQC → fastp
  → build human Bowtie2 index
  → Bowtie2 alignment to the human genome
  → unique + properly-paired filter → mito-% QC → keep chr1–22/chrX
  → Picard dedup → blacklist filter → MACS2 peaks
  → depth-normalized (RPGC) bigWigs
  → consensus peaks → featureCounts fragment matrix
```

**Key Features:**
- Quality assessment with FastQC; read trimming/filtering with fastp
- **Bowtie2 alignment** to the human genome
- Uniquely-mapped, properly-paired filtering; filtering-orphaned mates removed; **mitochondrial-% recorded** then non-primary chromosomes dropped
- Picard duplicate removal + fragment-level ENCODE blacklist filtering
- MACS2 peak calling (BAMPE)
- **Depth-normalized coverage**: RPGC (1×) bigWig tracks via deepTools `bamCoverage`
- **Consensus peaks** (fixed-width, SPM-ranked; majority-vote/IDR reproducibility) + **featureCounts** matrix

### 2. Quality Control Pipeline (`qc_all` target)

**Comprehensive QC Metrics** (run after the primary pipeline; consumes its `results/`):
- **Fragment Size Analysis** - Insert size distribution and nucleosomal patterns
- **TSS Enrichment** - Heatmap/profile **and a numeric per-sample enrichment score**
- **Signal Quality** - Fingerprint plots for signal-to-noise assessment
- **Sample Correlation** - Multi-sample correlation heatmap/scatter and PCA
- **GC Bias Assessment** - Sequence composition bias evaluation
- **FRiP Scores** - Fraction of Reads in Peaks
- **Library Complexity** - PCR bottlenecking assessment (NRF, PBC1, PBC2)
- **IDR Analysis** - Irreproducible Discovery Rate on relaxed peak calls, per replicate pair
- **Peak + annotation summary** - Peak counts/widths and reads in promoters vs enhancers
- **Interactive HTML QC report** (`atacseq_qc_report.html`) - one self-contained, theme-aware page for all QC **except FastQC** (alignment rate, mito %, duplication, blacklist, peaks/FRiP, TSS enrichment, NRF/PBC1/PBC2, fragment size/GC/correlation/PCA/fingerprint, reads-in-annotation, consensus), with ENCODE-threshold pass/warn/fail flags; numeric metrics render as interactive tables/bar charts and the deepTools QC (fragment size, GC bias, correlation/PCA, fingerprint, TSS profile + per-region heatmap) renders as **interactive SVG/canvas charts drawn client-side** — nothing is embedded as a static image. Mitochondrial % comes from the primary pipeline's `idxstats`
- **FastQC MultiQC report** (`multiqc_fastqc.html`) - MultiQC scoped to FastQC (raw-read quality) only

### 3. Differential Binding Analysis (`ATACseq_Dx.ipynb`)

An **R / Bioconductor** notebook (runs on the `atacseq-diffbind` container, `ir` kernel) that
takes the pipeline's consensus fragment matrix (`results/consensus/consensus_counts.txt`) and
performs the differential test end-to-end:

**Differential binding (DESeq2):**
- **NICD3 vs Ctrl** on the consensus peaks, split into **promoter vs distal** sets (overlap with
  `ref/promoter_chr1-22X.bed`) and tested separately
- **Paired design** `~ pair + condition` — each replicate index blocks pair-to-pair variation
- Normalized with DESeq2's **default** median-of-ratios size factors (`sf = NULL`)
- Per group: full results + significant subset (`padj < 0.05 & |log2FC| > 1`),
  MA/volcano plots, sample **PCA** (VST), and ChIPseeker nearest-gene
  annotation → `results/diff_region/`

**Positive-control genome-browser tracks (Gviz):**
- Signal at canonical Notch/NICD targets (HES1, HEY1, HEYL, NRARP, DTX1) across all samples, from
  the **RPGC** bigWigs, with all-transcript
  gene models from the GTF → `results/browser_tracks/` (PNG + PDF)

The executed notebook is generated from `workflow/scripts/build_diffbind_notebook.py`
(helpers in `workflow/scripts/diffbind_helpers.R`).

## Requirements

The pipeline requires the following dependencies:

- [Snakemake](https://snakemake.readthedocs.io/) ≥7.0.0
- [Conda](https://docs.conda.io/en/latest/) / [Mamba](https://github.com/mamba-org/mamba) (recommended)
- [Python](https://www.python.org/) ≥3.8
- UNIX-based system (Linux/MacOS)

### Software Dependencies
(automatically installed via conda environments):
- **FastQC** (quality control)
- **fastp** (read trimming)
- **Bowtie2** (alignment to the human genome)
- **SAMtools** (BAM processing, read splitting, idxstats)
- **Picard** (duplicate removal)
- **MACS2** (peak calling)
- **deepTools** (bigWig/bedGraph, QC and visualization)
- **bedtools** + **bc** (genomic interval operations; FRiP / complexity / annotation)
- **IDR** (reproducibility analysis)
- **featureCounts / Subread** (fragment quantification over the consensus set)

### Differential-analysis environment (R / Bioconductor)
`ATACseq_Dx.ipynb` is an **R** notebook (`ir` kernel), not Python. Run it on the
**`atacseq-diffbind`** container, which bundles R 4.6 / Bioconductor with everything it needs:
- **DESeq2**, **apeglm** — differential binding + log2FC shrinkage
- **GenomicRanges** / **IRanges**, **rtracklayer** — region handling, bigWig I/O
- **ChIPseeker**, **TxDb.Hsapiens.UCSC.hg38.knownGene**, **org.Hs.eg.db** — annotation
- **Gviz** — genome-browser tracks; **ggplot2** — MA / volcano / PCA

(edgeR / limma / DiffBind are also installed for alternative differential frameworks.)

### Reference Files (`ref/`)

Files fall into three groups: **shipped** with the repo (present after clone),
**downloaded** from public sources, and **generated** locally from the downloads.
Exact commands are in [Reference files: download & generate](#reference-files-download--generate).

```
ref/
├── config.yaml                          # pipeline configuration                 (shipped)
├── samples.csv                          # sample sheet: sample_id,type,group      (shipped; you edit)
│
│  ── downloaded (public sources; not in the repo) ──
├── hg38.fa                              # human genome FASTA (UCSC, chr-prefixed)
├── hg38.2bit                            # human genome 2bit (QC: GC bias)
├── gencode.v36.annotation.gtf           # GENCODE annotation (TSS, gene models)
├── FANTOM5_CAGE_peaks_hg38.bed.gz       # FANTOM5 CAGE peaks (optional; only to rebuild the CAGE set)
├── picard.jar                           # Picard (MarkDuplicates)
│
│  ── shipped, or generated from the downloads ──
├── hg38_blacklist_regions.bed           # ENCODE hg38 blacklist v2                (shipped)
├── hg38.fa.fai                          # faidx of hg38.fa                        (generated)
├── promoter_chr1-22X.bed                # Ensembl Regulatory Build promoters      (shipped; QC reads-in-annotation)
├── enhancer_chr1-22X.bed                # Ensembl Regulatory Build enhancers      (shipped; QC reads-in-annotation)
├── Promoter_uniqueTSS_hg38_{3000,5000}bp_chr1-22X.bed      # TSS±N, unique TSS          (shipped)
├── Promoter_MANEcanonical_hg38_{3000,5000}bp_chr1-22X.bed  # TSS±N, one per gene (MANE) (shipped)
├── Promoter_FANTOM5CAGE_hg38_{3000,5000}bp_chr1-22X.bed    # TSS±N, empirical CAGE      (shipped)
│      # per-transcript set (Promoter_UCSC_hg38_*bp) not shipped — `build_promoter_beds.py transcript`
├── BOWTIE2/                             # human Bowtie2 index                     (built by the pipeline)
│
│  ── shipped scripts ──
├── build_promoter_beds.py               # regenerates the promoter TSS BEDs (4 definitions)
├── build_diffbind_notebook.py           # regenerates ATACseq_Dx.ipynb
├── build_qc_report.py  tss_score.py  consensus_peaks.py  process_sam.py
├── blacklist-stats-script.py  downsample_tss_matrix.py
└── diffbind_helpers.R                   # R helpers for the differential-binding notebook
```

### Reference files: download & generate

Only the **downloaded** files are absent after a clone. Fetch them into `ref/`, then
build the small derived files. Everything else ships with the repo.

**1. Download** (run from the repo root):

```bash
cd ref

# Human genome and the human 2bit (UCSC goldenPath)
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz && gunzip hg38.fa.gz
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit

# GENCODE v36 gene annotation
curl -O http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gtf.gz
gunzip gencode.v36.annotation.gtf.gz

# FANTOM5 robust CAGE peaks (hg38) — OPTIONAL: only needed to rebuild the CAGE
# promoter set (the built Promoter_FANTOM5CAGE_*.bed already ships with the repo)
curl -L -o FANTOM5_CAGE_peaks_hg38.bed.gz \
  "https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/CAGE_peaks/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz"

# Picard (MarkDuplicates)
curl -L -o picard.jar https://github.com/broadinstitute/picard/releases/latest/download/picard.jar

cd ..
```

The ENCODE blacklist ships with the repo; to refresh it (Boyle Lab v2):

```bash
curl -L https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz \
  | gunzip > ref/hg38_blacklist_regions.bed
```

**2. Generate** (needs the downloads above; run from the repo root):

```bash
# faidx of the human genome (used to clamp promoter windows to chromosome ends)
samtools faidx ref/hg38.fa

# The promoter TSS BEDs (unique / MANE / CAGE) already ship with the repo. Rebuild
# them only after updating the GTF, or to also produce the per-transcript set:
python workflow/scripts/build_promoter_beds.py all   # modes: unique | mane | cage | transcript | all
# (needs ref/hg38.fa.fai + the GTF; the `cage` mode also needs the FANTOM5 download above)
# Kept chromosomes default to `keep_chroms` in config/config.yaml (currently chr1-22,X);
# override per run with e.g. --chroms chr1,chr2,chrX  or  --chroms all  (the filename
# scope token, e.g. chr1-22X, auto-adjusts; set it explicitly with --label).

# The human Bowtie2 index (ref/BOWTIE2/) is built automatically
# by the pipeline's build_genome_index rule from hg38.fa — no manual step.
```

The QC annotation BEDs (`promoter_chr1-22X.bed`, `enhancer_chr1-22X.bed`) ship with the
repo. They are **Ensembl Regulatory Build** features (`feature_type` Promoter / Enhancer,
`ENSR*` IDs), chr-prefixed and restricted to chr1–22,X. To rebuild from a newer release,
download the regulatory GFF from
`http://ftp.ensembl.org/pub/current_regulation/homo_sapiens/`, keep the Promoter /
Enhancer features, and emit `chrom start end ENSR_id`.

### Promoter / TSS definitions

`build_promoter_beds.py` writes promoter windows (TSS ±3kb and ±5kb, strand-aware,
chr1–22,X) under four TSS definitions so you can match the definition to the analysis.
All windows are symmetric width 2N; GTF-derived sets keep a 10-column layout
(`tx_id chrom start end score strand tx_id.ver transcript_type gene_name gene_id`),
CAGE is BED6.

| `mode` | TSS = | Windows (chr1–22,X) | Use when |
|---|---|---|---|
| `transcript` | every transcript's 5′ end (GENCODE v36) | 231,104 | you want all annotated alternative promoters (one per isoform; ChIPseeker-style) |
| `unique` | distinct (chrom, TSS, strand) tuples | 205,696 | as `transcript`, without isoform double-counting |
| `mane` | one canonical TSS per gene (MANE Select, else 5′-most) | 60,058 | one row per gene; ENCODE-comparable TSS-enrichment reference |
| `cage` | FANTOM5 robust CAGE peak position (empirical) | 209,146 | data-defined start sites, incl. promoters annotation misses |

The repo ships the **`unique`, `mane`, and `cage`** sets (each at ±3kb and ±5kb); the
per-transcript set is generate-on-demand (`build_promoter_beds.py transcript`).

Cross-check: **90.6%** of protein-coding MANE-canonical TSSs have a FANTOM5 CAGE peak
within 500 bp (strong annotation↔empirical agreement); non-coding genes only **24.9%**
(lncRNAs etc. are largely CAGE-silent). Note that each transcript has exactly one TSS,
whereas a gene can have several — the pipeline's deepTools TSS-enrichment QC and the
`ChIPseeker` peak annotation both use **transcript-level** TSS by default.

## Installation

```bash
# Clone the repository
git clone https://github.com/gynecoloji/snakemake_ATACseq.git
cd snakemake_ATACseq

# You need Snakemake + conda/mamba as the driver. The per-rule tool environments
# (workflow/envs/*.yaml) are created automatically on the first `--use-conda` run —
# you do not build them by hand.
mamba create -n atacseq -c conda-forge -c bioconda snakemake-minimal pandas
conda activate atacseq
```

> **No local install?** Skip all of the above and use the container image instead —
> see [Container Execution (Docker / Apptainer)](#container-execution-docker--apptainer).

## Configuration

All parameters live in `config/config.yaml`, which ships with working defaults and
an inline comment on each one. The [config schema](workflow/schemas/config.schema.yaml)
is the **single source of truth** for parameter types, defaults, and descriptions:
the workflow validates your config against it on every run (and fills in defaults
for anything you omit), and the
[Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=gynecoloji/snakemake_ATACseq)
renders it as a parameter table on the workflow page. See
[`config/README.md`](config/README.md) for the sample sheet and reference-data details.

At minimum, point the reference-file paths at the genomes/annotation you provide
(see [Data Preparation](#data-preparation)):

```yaml
samples_table: "config/samples.csv"                # sample sheet (sample_id, type, group)
human_fasta:   "ref/hg38.fa"                       # chr-prefixed UCSC hg38 (you provide)
gtf:           "ref/gencode.v36.annotation.gtf"    # GENCODE, chr-prefixed (for TSS-enrichment QC)
```

## Data Preparation

### Input Files
Please pay attention to the **common suffix** of the fastq related raw files (_R1/2_001.fastq.gz)
Accepted file format should look like: **sample_id**+**common suffix** (e.g. GSF4007-Control_1_S11_R1_001.fastq.gz)

Place paired-end FASTQ files in the `data/` directory following this naming convention:
```
data/{sample}_R1_001.fastq.gz
data/{sample}_R2_001.fastq.gz
```

### Sample Information

Create `config/samples.csv`:
```csv
sample_id,type,group
GSF4007-Control_1_S11,Control,group1
GSF4007-Control_2_S13,Control,group1
GSF4007-Control_3_S15,Control,group1
GSF4007-NICD3-V5_1_S12,NICD3,group2
GSF4007-NICD3-V5_2_S14,NICD3,group2
GSF4007-NICD3-V5_3_S16,NICD3,group2
```

The `group` column defines conditions / replicate sets. It drives **consensus-peak
reproducibility** (≥2-of-N majority vote for ≥3-replicate conditions, or IDR for
2-replicate conditions) and pairwise IDR in the QC pipeline.

## Running the Pipeline

The workflow is the standard-layout `workflow/Snakefile`. A single run builds the
primary stage **and** the QC report in dependency order (unified DAG); use the
`atacseq_all` / `qc_all` targets to run just one stage. Snakemake auto-discovers
`workflow/Snakefile` when run from the repo root, so `-s workflow/Snakefile` is
optional but shown here for clarity.

### Dry Run

To check the workflow without executing any commands:
```bash
# Check the whole workflow (primary → QC)
snakemake -s workflow/Snakefile -n

# Or check a single stage
snakemake -s workflow/Snakefile -n atacseq_all
snakemake -s workflow/Snakefile -n qc_all
```

### Local Execution

```bash
# Run everything (primary stage → QC report) in one dependency-ordered DAG
snakemake -s workflow/Snakefile --use-conda --cores 20

# Or run a single stage
snakemake -s workflow/Snakefile --use-conda --cores 20 atacseq_all   # primary only
snakemake -s workflow/Snakefile --use-conda --cores 20 qc_all        # QC only (after primary)
```

### Differential Analysis
`ATACseq_Dx.ipynb` is an **R** notebook (`ir` kernel); run it in the **`atacseq-diffbind`**
environment. Interactively:
```bash
jupyter notebook ATACseq_Dx.ipynb
```
Or execute it headless (writes `results/diff_region/` + `results/browser_tracks/`):
```bash
jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.kernel_name=ir ATACseq_Dx.ipynb
```
No local R/Bioconductor install? Run it in the **`atacseq-diffbind`** container (pull it
from [Container Execution](#container-execution-docker--apptainer)); Apptainer auto-mounts
the current directory:
```bash
apptainer exec atacseq-diffbind.sif \
    jupyter nbconvert --to notebook --execute --inplace \
    --ExecutePreprocessor.kernel_name=ir ATACseq_Dx.ipynb
```

### Cluster Execution

For execution on a SLURM cluster: (Not tested)
```bash
snakemake -s workflow/Snakefile --use-conda \
  --cluster "sbatch -p {params.partition} -c {threads} -t {params.time}" \
  --jobs 100
```

### Container Execution (Docker / Apptainer)

Two prebuilt images cover the whole workflow — you install nothing except Docker or
Apptainer:

| Image | Contents | Used for |
|---|---|---|
| **`gynecoloji/atacseq`** | Snakemake + all five per-rule conda envs | primary pipeline + QC pipeline |
| **`gynecoloji/atacseq-diffbind`** | R / Bioconductor (DESeq2, ChIPseeker, Gviz) + `ir` Jupyter kernel | `ATACseq_Dx.ipynb` differential analysis |

**Download** — pull with Docker, or convert to a local `.sif` once for Apptainer /
Singularity (HPC):

```bash
# Docker
docker pull gynecoloji/atacseq:latest
docker pull gynecoloji/atacseq-diffbind:latest

# Apptainer / Singularity  (writes ./atacseq-*.sif in the current directory)
apptainer pull atacseq.sif  docker://gynecoloji/atacseq:latest
apptainer pull atacseq-diffbind.sif docker://gynecoloji/atacseq-diffbind:latest
```

Genomes/FASTQs are **not** baked into the images; you mount your project directory at run
time (see [`DOCKER.md`](DOCKER.md) for the exact `ref/` and `data/` files the container
expects). A single run builds the primary stage then QC (unified DAG); the differential
notebook runs in the diffbind image (see [Differential Analysis](#differential-analysis)).

The image's entrypoint is
`snakemake --use-conda --conda-frontend mamba --conda-prefix /opt/wf-conda`, so anything
you pass after the image name goes straight to `snakemake` (e.g. `-s workflow/Snakefile
--cores N`, an optional `atacseq_all`/`qc_all` target, or `-n` for a dry run).

#### Docker

```bash
# Pull the published image (or run `docker compose build` to build it locally)
docker pull gynecoloji/atacseq:latest

# Run from your project directory (which holds workflow/, config/, ref/, data/):
# Everything (primary stage → QC report) in one dependency-ordered run:
docker run --rm -v "$(pwd)":/workflow -e HOME=/tmp --user "$(id -u):$(id -g)" \
    gynecoloji/atacseq:latest -s workflow/Snakefile --cores 16

# Or just one stage: append the atacseq_all or qc_all target
docker run --rm -v "$(pwd)":/workflow -e HOME=/tmp --user "$(id -u):$(id -g)" \
    gynecoloji/atacseq:latest -s workflow/Snakefile --cores 16 qc_all

# dry run: append  -n
```

Convenience wrappers `docker compose` and `./run_pipeline.sh` are also provided:

```bash
./run_pipeline.sh --cores 16                 # everything (primary → QC)
docker compose run --rm atacseq --cores 16 qc_all
```

#### Apptainer / Singularity (HPC)

On clusters without Docker, convert the image to a SIF once and run it with Apptainer.
Apptainer auto-mounts `$HOME`, `/tmp`, and the current directory, and runs as you (no
`--user` needed):

```bash
# One-time: build a local .sif from the Docker Hub image
apptainer pull atacseq.sif docker://gynecoloji/atacseq:latest

# Run from your project directory (everything: primary stage → QC report):
apptainer exec atacseq.sif \
    snakemake --use-conda --conda-frontend mamba --conda-prefix /opt/wf-conda \
    -s workflow/Snakefile --cores 16

# Or just one stage: append the atacseq_all or qc_all target
apptainer exec atacseq.sif \
    snakemake --use-conda --conda-frontend mamba --conda-prefix /opt/wf-conda \
    -s workflow/Snakefile --cores 16 qc_all
```

Notes:
- **References/data outside the project dir:** if `ref/` genomes live elsewhere (e.g. on
  scratch), bind them in — `--bind /scratch/genomes:/scratch/genomes` — and point
  `config/config.yaml` at the bound paths.
- **Pre-built envs:** the five conda environments are baked at `/opt/wf-conda`
  (read-only in the SIF) and reused via `--conda-prefix`. If Apptainer reports a
  read-only error writing there, add `--writable-tmpfs` to the `apptainer exec` command.
- `apptainer run atacseq.sif -s workflow/Snakefile --cores 16` also works — it
  invokes the same entrypoint.

## Deploying with snakedeploy

This repository follows the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/)
standardized structure (`workflow/Snakefile`, `config/`, `workflow/rules|scripts|envs/`,
and `.snakemake-workflow-catalog.yml`), so it can be deployed into another project
without cloning it by hand:

```bash
pip install snakedeploy
# In an empty target project directory:
snakedeploy deploy-workflow https://github.com/gynecoloji/snakemake_ATACseq . --tag main
```

This writes `workflow/Snakefile` (which `module`-imports this workflow) and a
`config/` copy for you to edit. You can also import selected rules into your own
`Snakefile` with Snakemake's module system:

```python
module atacseq:
    snakefile:
        github("gynecoloji/snakemake_ATACseq", path="workflow/Snakefile", tag="main")
    config:
        config

use rule * from atacseq
```

Then supply your own `config/config.yaml`, `config/samples.csv`, and `ref/` reference
data (see [Configuration](#configuration)) and run with `snakemake --use-conda`.

## Pipeline Details

### 1. Quality Control and Preprocessing

- **FastQC** - Quality assessment of raw reads
- **Fastp** - Adapter trimming and quality filtering with the following parameters:
  - Minimum read length: 30bp
  - Adapter handling: **auto-detects** adapters for paired-end reads by default
    (`--detect_adapter_for_pe`); set `adapter_r1`/`adapter_r2` in `config/config.yaml`
    to override with explicit sequences
  - Polyg tail trimming
  - Quality trimming: sliding window of 4 with mean quality 20

### 2. Genome Alignment

- **Genome index** (`build_genome_index`) - The human genome is optionally subset to
  `align_chroms`, and a Bowtie2 index is built at `bowtie2_index` (`ref/BOWTIE2/genome`).
- **Bowtie2** (`bowtie2_align`) - Alignment to the human genome:
  - `-X 3000 -I 0 --no-discordant --no-mixed`
- **Read filtering** (`samtools_sort_filter_index`):
  - Properly-paired (`-f 2 -F 2316`), uniquely-mapped (`grep -v XS:i:`),
    filtering-orphaned mates removed (`process_sam.py`)
- **Mitochondrial-% QC** - `samtools idxstats` on the BAM (incl. chrM) is recorded,
  then reads are restricted to `keep_chroms` (chr1–22, chrX; chrM/chrY/non-primary dropped).

### 3. Post-processing

- **Picard MarkDuplicates** - PCR duplicate removal with `REMOVE_DUPLICATES=true`
- **Fragment-level blacklist filtering** - ENCODE blacklist region removal using bedtools with proper paired-end handling

### 4. Peak Calling

- **MACS2** - Peak calling in paired-end mode:
  - Format: BAMPE, genome `hs`, `--nomodel`, q-value cutoff 0.05

### 5. Coverage Tracks and Consensus Peaks

- **Coverage tracks** (`create_bigwig`) - deepTools `bamCoverage` emits a depth-normalized
  RPGC (1×, "reads per genomic content") bigWig per sample for cross-sample comparison.
- **Consensus peaks** (`consensus_peaks`) - Corces-2018 fixed-width (`consensus_window`),
  score-per-million (SPM) iterative-overlap peaks. Per-condition reproducibility is a
  ≥2-of-N majority vote (≥3 replicates) or IDR (2 replicates, via `relaxed_peaks` +
  `reproducible_idr`), unioned across conditions, with blacklist / chrM / chrY / non-primary
  excluded.
- **Fragment counting** (`count_fragments_consensus`) - `featureCounts` (paired-end) over the
  consensus set → `results/consensus/consensus_counts.txt` (a regions × samples matrix).

### 6. Comprehensive QC Metrics (`qc_all` target)

- **Fragment Size Analysis** - Nucleosomal pattern assessment using `bamPEFragmentSize`
- **TSS Enrichment** - Heatmap/profile plus a **numeric enrichment score** per sample
- **Fingerprint Analysis** - Signal-to-noise using `plotFingerprint`
- **Sample Correlation** - `multiBamSummary` → correlation heatmap/scatter and PCA
- **GC Bias Assessment** - Sequence composition bias (`computeGCBias`)
- **Library Complexity**: **NRF** (Nd/Total), **PBC1** (N1/Nd), **PBC2** (N1/N2)
- **FRiP Score** - Fraction of reads in peaks
- **Peak + annotation summary** - Peak counts/widths; reads in promoters vs enhancers
- **IDR Analysis** - IDR on relaxed peak calls, per replicate pair
- **Interactive QC report** - `atacseq_qc_report.html` (self-contained; all QC except FastQC) plus `multiqc_fastqc.html` (FastQC only)
  (mitochondrial % comes from the primary pipeline's `idxstats`)

### 7. Differential Binding Analysis (`ATACseq_Dx.ipynb`)

An R / Bioconductor notebook (details in [Component 3](#3-differential-binding-analysis-atacseq_dxipynb)):
- **DESeq2** differential binding **NICD3 vs Ctrl** on the consensus matrix — **promoter vs distal**,
  **paired** design, under **default (median-of-ratios) normalization**
- Sample **PCA**, MA/volcano plots, ChIPseeker nearest-gene annotation → `results/diff_region/`
- **Gviz** genome-browser tracks at Notch target loci (RPGC bigWigs) →
  `results/browser_tracks/`

## Output Files

The pipeline generates the following output directories:

```
results/
├── fastqc/                 # FastQC reports
├── fastp/                  # Trimmed reads and reports
├── aligned/                # Bowtie2 log ({sample}.bowtie2.log); SAM is temporary
├── filtered/               # Filtered BAM + {sample}.idxstats.txt (mito QC) + summaries
├── dedup/                  # Deduplicated BAM + Picard metrics
├── blacklist_filtered/     # Analysis-ready BAM ({sample}.nobl.bam)
├── peaks/                  # MACS2 peak calls (*_peaks.narrowPeak)
├── bigwig/                 # Depth-normalized (RPGC) bigWigs ({sample}.bw)
├── consensus/
│   ├── consensus_peaks.bed # Fixed-width, non-overlapping consensus set
│   ├── consensus_peaks.saf # featureCounts input
│   └── consensus_counts.txt # Fragment count matrix (regions × samples)
├── peaks_relaxed/, qc_relaxed_peaks/  # Relaxed MACS2 calls (IDR input; 2-rep conditions / QC)
├── deeptools/              # Fragment size (+ raw lengths), fingerprint, correlation/PCA (+ matrix),
│                           #   GC bias, TSS matrix/plots + downsampled-heatmap JSON
├── bedgraph/               # RPGC bedGraphs (QC)
├── FRiP/                   # FRiP scores (*.frip.txt)
├── library_complexity/     # NRF / PBC1 / PBC2 (*_complexity.txt)
├── idr/                    # IDR peak calls between replicate pairs
├── peak_annotation/        # Reads in promoters / enhancers
├── qc/
│   ├── atacseq_qc_report.html      # Interactive QC report (all QC except FastQC)
│   ├── multiqc_fastqc.html         # FastQC-only MultiQC
│   ├── tss_enrichment_scores.tsv, peak_summary.tsv
│   └── blacklist_filtering_stats.txt
# ── Differential-analysis notebook (ATACseq_Dx.ipynb) outputs ──
├── diff_region/            # DESeq2 DB (promoter/distal), PCA, MA/volcano, annotation
└── browser_tracks/         # Gviz tracks at Notch loci (RPGC bigWigs, PNG + PDF)
```

### Directory Structure
```
ATAC-seq-Pipeline/                     # Snakemake Workflow Catalog layout
├── config/
│   ├── config.yaml         # Workflow parameters
│   ├── samples.csv         # Sample sheet (sample_id, type, group)
│   └── README.md           # Configuration reference
├── workflow/
│   ├── Snakefile           # Entry point (unified DAG; targets: atacseq_all, qc_all)
│   ├── rules/              # common.smk, atacseq.smk, qc.smk
│   ├── scripts/            # Python / R scripts used by the rules (+ helpers)
│   └── envs/               # Per-rule conda environment files
├── .snakemake-workflow-catalog.yml    # Catalog metadata (enables snakedeploy)
├── data/                   # Raw FASTQ files (you provide)
├── ref/                    # Reference genomes/annotations/BEDs (you download)
│   └── BOWTIE2/            # Human Bowtie2 index (built by the pipeline)
├── create_envs.smk         # Build-time helper (pre-bakes the conda envs into the image)
├── Dockerfile, docker-compose.yml, run_pipeline.sh, DOCKER.md
├── results/                # All pipeline outputs (detailed above)
│   └── tmp/                # Temporary processing files
└── logs/                   # Per-rule logs (see below)
```

## Troubleshooting

### Common Issues

1. **Low alignment rate** 
   - Check reference genome compatibility and version
   - Verify adapter trimming parameters in fastp step
   - Check for sample contamination or incorrect library prep
   - Examine FastQC reports for quality issues

2. **High duplication rate**
   - Indicates low library complexity issue
   - Consider optimizing chromatin extraction or Tn5 tagmentation conditions
   - May require increased sequencing depth for better coverage
   - Check PBC1 and PBC2 metrics for severity assessment

3. **Poor TSS enrichment**
   - Issues with sample quality or chromatin accessibility
   - Check protocol for nuclei isolation and tagmentation efficiency
   - Verify blacklist filtering is not removing genuine signal
   - Consider cell type-specific accessibility patterns

4. **High mitochondrial content**
   - Common in certain cell types but may indicate technical issues
   - Consider optimizing nuclear isolation protocols
   - May require additional sequencing depth to compensate for lost reads
   - Check for cytoplasmic contamination during sample prep

5. **Failed IDR analysis**
   - Low reproducibility between biological replicates
   - Check experimental conditions for consistency across replicates
   - May indicate technical issues with sample preparation or processing
   - Consider peak calling parameters if IDR consistently fails

6. **Fragment size distribution issues**
   - Should show nucleosomal ladder pattern (147bp, 294bp, etc.)
   - Flat distribution may indicate DNA degradation or poor tagmentation
   - Very short fragments may indicate over-tagmentation

### Performance Optimization
- Adjust thread counts based on available resources
- Use SSD storage for temporary files when possible
- Monitor memory usage during peak calling with large datasets
- Consider using `--cluster` mode for large-scale analyses

### Log Files

Comprehensive log files for each step are stored in the `logs/` directory:
```
logs/
# Primary pipeline (atacseq_all target)
├── fastqc/, fastp/, bowtie2/, samtools/, dedup/, blacklist_filter/, macs2/
├── build_genome_index/
├── bigwig/
├── relaxed_peaks/, reproducible/, consensus/, consensus_counts/
├── blacklist_stats/
# QC pipeline (qc_all target)
├── deeptools_bedgraph/, deeptools_fragmentsize/, deeptools_plotfingerprint/
├── deeptools_correlation/, deeptools_gc_bias/, deeptools_tss/, tss_enrichment_score/
├── FRiP/, qc_relaxed_peaks/, idr/, library_complexity/
└── reads_in_annotations/, peak_summary/, multiqc_fastqc/
```

## Citing the underlying tools

If you use this pipeline, please also cite the relevant tools:
- Snakemake: Köster & Rahmann (2012)
- MACS2: Zhang et al. (2008)
- deepTools: Ramírez et al. (2016)
- Bowtie2: Langmead & Salzberg (2012)
- IDR: Li et al. (2011); featureCounts: Liao et al. (2014); Consensus peaks: Corces et al. (2018)

## Contact

For questions or issues with this pipeline, please refer to the individual tool documentation or submit an issue to the repository.

---

**Note**: This pipeline is optimized for human genome analysis (hg38) but can be adapted for other organisms by updating reference files and parameters.