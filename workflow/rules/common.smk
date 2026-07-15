# Shared setup for the primary (atacseq.smk) and QC (qc.smk) rule files:
# imports, config-derived values, output-directory constants, sample/group
# tables and helper functions. Included first by workflow/Snakefile, so every
# name defined here is visible to the rules in the other included files.

import pandas as pd
import re
import os
from itertools import combinations
from snakemake.utils import validate

# ── Config validation ───────────────────────────────────────────────────
# Validate against workflow/schemas/config.schema.yaml (also fills in defaults
# for any omitted parameters). The catalog renders the parameter table from the
# same schema. Path is relative to this file (workflow/rules/).
validate(config, "../schemas/config.schema.yaml")

# ── Samples ─────────────────────────────────────────────────────────────
samples_df = pd.read_csv(config["samples_table"])
SAMPLES = samples_df["sample_id"].tolist()

# ── Output directories (all relative to the working dir) ────────────────
RESULT_DIR             = "results"
FASTQC_DIR             = f"{RESULT_DIR}/fastqc"
FASTP_DIR              = f"{RESULT_DIR}/fastp"
ALIGN_DIR              = f"{RESULT_DIR}/aligned"
TMP_DIR                = f"{RESULT_DIR}/tmp"
FILTERED_DIR           = f"{RESULT_DIR}/filtered"
DEDUP_DIR              = f"{RESULT_DIR}/dedup"
BLACKLIST_FILTERED_DIR = f"{RESULT_DIR}/blacklist_filtered"
PEAKS_DIR              = f"{RESULT_DIR}/peaks"
QC_DIR                 = f"{RESULT_DIR}/qc"

# ── Coverage (depth-normalized bigWig) directory ────────────────────────
BIGWIG_DIR         = f"{RESULT_DIR}/bigwig"

# ── Consensus (Module B) directories ────────────────────────────────────
RELAXED_PEAKS_DIR  = f"{RESULT_DIR}/peaks_relaxed"
CONSENSUS_DIR      = f"{RESULT_DIR}/consensus"

# ── QC-pipeline directories (aliases + QC-only outputs) ─────────────────
RMD_BAM_DIR    = BLACKLIST_FILTERED_DIR   # QC alias: results/blacklist_filtered
PEAK_DIR       = PEAKS_DIR                 # QC alias: results/peaks
BEDGRAPH_DIR   = f"{RESULT_DIR}/bedgraph"
DEEPTOOLS_DIR  = f"{RESULT_DIR}/deeptools"
FRIP_DIR       = f"{RESULT_DIR}/FRiP"
IDR_DIR        = f"{RESULT_DIR}/idr"
RELAXED_DIR    = f"{RESULT_DIR}/qc_relaxed_peaks"
COMPLEXITY_DIR = f"{RESULT_DIR}/library_complexity"
ANNOT_DIR      = f"{RESULT_DIR}/peak_annotation"

# ── Reference data / config ─────────────────────────────────────────────
GENOME_2BIT  = os.path.join("ref", "hg38.2bit")   # QC: computeGCBias --genome
GTF_FILE     = config["gtf"]
PROMOTER_BED = config["promoter_bed"]
ENHANCER_BED = config["enhancer_bed"]
EGS          = config["effective_genome_size"]
PEAK_TYPES   = config.get("peak_types", ["narrowPeak"])

# ── Conditions (groups) and per-condition reproducibility method ────────
GROUPS = samples_df.groupby("group")["sample_id"].apply(list).to_dict()

def _repro_method(members):
    n = len(members)
    if n >= 3:
        return "majority"
    if n == 2:
        return "idr"
    return "single"

GROUP_METHOD  = {g: _repro_method(m) for g, m in GROUPS.items()}
IDR_GROUPS    = [g for g in GROUPS if GROUP_METHOD[g] == "idr"]
NONIDR_GROUPS = [g for g in GROUPS if GROUP_METHOD[g] != "idr"]
IDR_SAMPLES   = [s for g in IDR_GROUPS for s in GROUPS[g]]

# IDR pairs (QC): all within-condition replicate pairs
IDR_PAIRS = []
for group, members in GROUPS.items():
    for a, b in combinations(members, 2):
        IDR_PAIRS.append((group, a, b))

def _alt(names):
    """Regex alternation for wildcard_constraints; matches nothing if empty."""
    return "|".join(re.escape(n) for n in names) if names else "a^"

# fastp adapter handling: AUTO-DETECT adapters for paired-end reads by default
# (--detect_adapter_for_pe). If adapter sequences are provided in config
# (adapter_r1 / adapter_r2, non-empty), pass them explicitly instead — that
# OVERRIDES auto-detection. Leave them unset/empty to auto-detect.
def _fastp_adapter_args():
    r1 = str(config.get("adapter_r1") or "").strip()
    r2 = str(config.get("adapter_r2") or "").strip()
    if r1:
        args = f"--adapter_sequence {r1}"
        if r2:
            args += f" --adapter_sequence_r2 {r2}"
        return args
    return "--detect_adapter_for_pe"

FASTP_ADAPTER_ARGS = _fastp_adapter_args()

def _group_relaxed_inputs(wildcards):
    return [f"{RELAXED_PEAKS_DIR}/{s}_relaxed.narrowPeak" for s in GROUPS[wildcards.group]]
