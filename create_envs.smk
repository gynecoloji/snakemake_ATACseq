# Build-time helper: pre-create the per-rule conda envs (envs/*.yaml) so they are
# baked into the Docker image and reused at runtime via --conda-prefix.
#
# Used ONLY by the Dockerfile:
#   snakemake -s create_envs.smk --use-conda --conda-create-envs-only \
#       --conda-frontend mamba --conda-prefix /opt/wf-conda --cores 1
#
# It has no external inputs, so the DAG resolves without any genomes/FASTQs.
# Snakemake keys each conda env by the CONTENT of its workflow/envs/*.yaml file,
# so the envs built here are reused by workflow/Snakefile at runtime (same env
# files, same --conda-prefix). If content differs, Snakemake just rebuilds that
# one env at runtime — no hard failure.

ENVS = ["snakemake", "deeptools", "macs2", "idr", "bedtools"]

rule all:
    input:
        expand("build/conda_env_{env}.ready", env=ENVS)

rule create_env:
    output:
        "build/conda_env_{env}.ready"
    conda:
        "workflow/envs/{env}.yaml"
    shell:
        "touch {output}"
