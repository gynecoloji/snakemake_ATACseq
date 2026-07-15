#!/usr/bin/env bash
# Build (once) and run the ATAC-seq Snakemake workflow in Docker.
#
# The workflow is the standard-layout workflow/Snakefile: a single run builds the
# primary pipeline AND the QC report (unified DAG). Pass extra snakemake args
# (cores, a target, -n, ...) through to the container.
#
# Usage:
#   ./run_pipeline.sh                          # everything (primary → QC), 4 cores
#   ./run_pipeline.sh --cores 16               # everything, 16 cores
#   ./run_pipeline.sh --cores 16 atacseq_all   # primary pipeline only
#   ./run_pipeline.sh --cores 16 qc_all        # QC only (after primary)
#   ./run_pipeline.sh -n                        # dry run: check the DAG first
#
# The current directory (code + ref/ genomes + data/ FASTQs) is mounted at
# /workflow; results/ are written back here. The pre-built conda envs live in the
# image at /opt/wf-conda and are reused via --conda-prefix. The image ENTRYPOINT
# runs `snakemake -s workflow/Snakefile --use-conda ...`.
set -euo pipefail

IMAGE="atacseq:latest"

if [ "$#" -eq 0 ]; then set -- --cores 4; fi   # default snakemake args

if ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
    echo ">> Building $IMAGE (first time only; pre-builds the conda envs, ~15-30 min)..."
    docker build -t "$IMAGE" .
fi

echo ">> snakemake -s workflow/Snakefile $*"
docker run --rm \
    -v "$(pwd)":/workflow \
    -e HOME=/tmp \
    --user "$(id -u):$(id -g)" \
    "$IMAGE" -s workflow/Snakefile "$@"
