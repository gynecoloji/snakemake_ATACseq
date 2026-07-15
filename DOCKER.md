# Running the ATAC-seq workflow with Docker

The pipeline uses **one conda environment per rule** (`workflow/envs/*.yaml`) because its
tools require incompatible Python versions (`idr`=3.6, `macs2`=3.7,
`snakemake`/`deeptools`=3.12). The image therefore ships **Snakemake + the 5
pre-built conda envs** and runs Snakemake with `--use-conda`.

Large reference genomes and FASTQs are **not** baked into the image — you mount
your project directory at run time.

## 1. What you need on the host (in `ref/` and `data/`)

The container reads these from the mounted project directory:

| Path | What |
|---|---|
| `data/{sample}_R1_001.fastq.gz`, `_R2_001.fastq.gz` | your paired-end reads |
| `ref/hg38.fa` | human genome FASTA (chr-prefixed UCSC) |
| `ref/hg38_blacklist_regions.bed` | ENCODE blacklist |
| `ref/picard.jar` | Picard (used by dedup) |
| `ref/gencode.v36.annotation.gtf`, `ref/hg38.2bit` | QC (TSS, GC bias) |
| `ref/promoter_chr1-22X.bed`, `ref/enhancer_chr1-22X.bed` | QC (reads-in-annotation) |
| `config/config.yaml`, `config/samples.csv` | config + sample sheet (tracked in the repo) |

The Bowtie2 index is built by the pipeline itself (`build_genome_index`).

## 2. Build the image (once)

```bash
docker compose build
# or:  docker build -t atacseq:latest .
```

This pre-builds the 5 conda envs into the image (a few GB; ~15–30 min the first
time). For a reproducible image, pin the base tag in the `Dockerfile`
(`FROM condaforge/miniforge3:<version>`).

## 3. Run

A single run builds the primary pipeline **and** the QC report (unified DAG).
Using the helper script (recommended):

```bash
./run_pipeline.sh -n                        # dry run: check the DAG first
./run_pipeline.sh --cores 16                # everything (primary → QC)
./run_pipeline.sh --cores 16 atacseq_all    # primary pipeline only
./run_pipeline.sh --cores 16 qc_all         # QC only (after primary)
```

Or with docker compose (the image entrypoint sets `-s workflow/Snakefile`; run
from the project root so Snakemake finds `workflow/Snakefile`):

```bash
docker compose run --rm atacseq -n
docker compose run --rm atacseq --cores 16
docker compose run --rm atacseq --cores 16 atacseq_all
docker compose run --rm atacseq --cores 16 qc_all
```

Or a raw `docker run` (mount the project; reuse the baked envs):

```bash
docker run --rm -v "$(pwd)":/workflow -e HOME=/tmp --user "$(id -u):$(id -g)" \
    atacseq:latest -s workflow/Snakefile --cores 16
```

Everything after the image name is passed straight to `snakemake` (the image's
entrypoint already sets `--use-conda --conda-frontend mamba --conda-prefix
/opt/wf-conda`).

**Targets:** the default target runs primary → QC in dependency order. Use the
`atacseq_all` target for just the primary pipeline (alignment → peaks → bigWigs →
consensus) and `qc_all` for just the QC stage (it consumes the primary pipeline's
`results/`).

## 4. Notes & troubleshooting

- **First run builds the Bowtie2 index** (`ref/BOWTIE2/…`) from `hg38.fa` — a
  large one-time step. It's cached for later runs.
- **Outputs ownership:** the run script / compose run as your host UID/GID
  (`--user`) so `results/` isn't root-owned. For compose, export `DOCKER_UID`/
  `DOCKER_GID` if the defaults (1000:1000) aren't you.
- **`defaults` channel ToS:** the env YAMLs list the Anaconda `defaults` channel.
  The Dockerfile best-effort-accepts its ToS; if an env solve still fails on
  `defaults`, either accept it (`conda tos accept …`) or drop `- defaults` from
  the affected `workflow/envs/*.yaml`.
- **`bc` for FRiP/complexity:** those QC rules use `workflow/envs/bedtools.yaml`, which
  must contain `bc` (and `samtools`, `bedtools`).
- **Cores:** pass `--cores N` to match the host; add `--resources mem_mb=…` if you
  cap memory. The `bowtie2-build` and alignments are the heavy steps.
