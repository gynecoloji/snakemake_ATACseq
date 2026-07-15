# ATAC-seq (depth-normalized) + consensus + QC — Snakemake workflow container.
#
# The workflow runs one conda env PER RULE (envs/*.yaml) because its tools need
# incompatible Python versions (idr=3.6, macs2=3.7, snakemake/deeptools=3.12).
# So this image ships Snakemake + the 5 pre-built conda envs and runs --use-conda.
#
# Large genomes/FASTQs are NOT baked in — mount your project directory at runtime
# (see docker-compose.yml / run_pipeline.sh / DOCKER.md).

# For a fully reproducible build, pin to a specific tag, e.g.
#   FROM condaforge/miniforge3:24.11.3-2
FROM condaforge/miniforge3:latest

LABEL org.opencontainers.image.title="atacseq-snakemake"
LABEL org.opencontainers.image.description="ATAC-seq (depth-normalized) + consensus + QC (Snakemake, --use-conda)"

ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    WF_CONDA_PREFIX=/opt/wf-conda

# Minimal system deps (git helps Snakemake's conda handling; procps for subprocess mgmt)
RUN apt-get update && apt-get install -y --no-install-recommends \
        git procps ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# Use FLEXIBLE channel priority: the env files are fully-pinned exports whose exact
# builds come from a mix of conda-forge/bioconda/defaults, so strict priority would
# wrongly exclude pins that live in a lower-priority channel. Also accept the
# Anaconda "defaults" ToS so the non-interactive solve doesn't stall.
RUN conda config --system --set channel_priority flexible && \
    ( conda tos accept --override-channels \
        --channel https://repo.anaconda.com/pkgs/main \
        --channel https://repo.anaconda.com/pkgs/r 2>/dev/null || true )

# Snakemake driver in its OWN env: the miniforge base pins Python 3.13, but
# snakemake-minimal 9.3.2 needs Python <3.13. pandas is required because the
# Snakefiles import it at parse time. Prepend the driver env to PATH.
RUN mamba create -y -n driver -c conda-forge -c bioconda \
        python=3.12 snakemake-minimal=9.3.2 pandas && \
    mamba clean -afy
ENV PATH=/opt/conda/envs/driver/bin:$PATH

WORKDIR /workflow

# Env specs first, then pre-build the 5 per-rule conda envs INTO the image. Doing
# this BEFORE copying the workflow code means later edits to the Snakefiles/scripts
# don't invalidate the (slow) conda-env layer.
COPY workflow/envs/ ./workflow/envs/
COPY create_envs.smk ./
RUN snakemake -s create_envs.smk --use-conda --conda-create-envs-only \
        --conda-frontend mamba --conda-prefix "${WF_CONDA_PREFIX}" --cores 1 && \
    mamba clean -afy && \
    rm -rf build .snakemake

# Workflow code + config (genomes/data are mounted at runtime, not baked)
COPY workflow/ ./workflow/
COPY config/ ./config/
COPY tests/ ./tests/

# ENTRYPOINT fixes the conda settings; pass the target/cores at `docker run`.
ENTRYPOINT ["snakemake", "--use-conda", "--conda-frontend", "mamba", "--conda-prefix", "/opt/wf-conda"]
CMD ["-s", "workflow/Snakefile", "--cores", "4"]
