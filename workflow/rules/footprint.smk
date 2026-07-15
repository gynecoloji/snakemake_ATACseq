# Optional TOBIAS differential TF footprinting stage. Reuses the primary
# pipeline's blacklist-filtered BAMs and the Module B consensus peaks; nothing
# here runs unless you request the footprint_all target:
#
#   snakemake --use-conda --cores N footprint_all
#
# Conditions are the distinct values of the sample sheet's `type` column
# (e.g. Control, NICD3). Per condition, the replicate BAMs are pooled, Tn5
# insertion bias is modelled and removed (ATACorrect — the rigorous alternative
# to a fixed +4/-5 shift), footprint scores are computed (ScoreBigwig), and
# BINDetect scans the JASPAR motifs for differential TF binding between
# conditions. Requires a JASPAR-format motif file at config['jaspar_motifs'].

CONDITION_NAMES = list(CONDITIONS.keys())


wildcard_constraints:
    cond = _alt(CONDITION_NAMES)


# Aggregate target for the footprinting stage (opt-in; NOT part of `rule all`).
rule footprint_all:
    input:
        expand(f"{FOOTPRINT_DIR}/{{cond}}.merged.bam", cond=CONDITION_NAMES),
        expand(f"{FOOTPRINT_DIR}/atacorrect/{{cond}}_corrected.bw", cond=CONDITION_NAMES),
        expand(f"{FOOTPRINT_DIR}/{{cond}}_footprints.bw", cond=CONDITION_NAMES),
        f"{FOOTPRINT_DIR}/bindetect/bindetect_results.txt"


def _condition_bams(wildcards):
    """Blacklist-filtered BAMs for every sample of this condition (`type`)."""
    return [f"{BLACKLIST_FILTERED_DIR}/{s}.nobl.bam" for s in CONDITIONS[wildcards.cond]]


# Pool the replicate BAMs of a condition into one signal for footprinting.
rule footprint_merge_condition:
    input:
        bams = _condition_bams
    output:
        bam = f"{FOOTPRINT_DIR}/{{cond}}.merged.bam",
        bai = f"{FOOTPRINT_DIR}/{{cond}}.merged.bam.bai"
    threads: 8
    conda:
        "../envs/tobias.yaml"
    log:
        "logs/footprint_merge/{cond}.log"
    shell:
        """
        mkdir -p {FOOTPRINT_DIR} logs/footprint_merge
        samtools merge -f -@ {threads} {output.bam} {input.bams} 2> {log}
        samtools index -@ {threads} {output.bam} 2>> {log}
        """


# TOBIAS ATACorrect: model + remove Tn5 insertion bias over the accessible
# (consensus) regions, producing a bias-corrected signal bigWig per condition.
rule tobias_atacorrect:
    input:
        bam       = f"{FOOTPRINT_DIR}/{{cond}}.merged.bam",
        bai       = f"{FOOTPRINT_DIR}/{{cond}}.merged.bam.bai",
        genome    = config["human_fasta"],
        peaks     = f"{CONSENSUS_DIR}/consensus_peaks.bed",
        blacklist = config["blacklist"]
    output:
        corrected = f"{FOOTPRINT_DIR}/atacorrect/{{cond}}_corrected.bw"
    params:
        outdir = f"{FOOTPRINT_DIR}/atacorrect"
    threads: 16
    conda:
        "../envs/tobias.yaml"
    log:
        "logs/tobias_atacorrect/{cond}.log"
    shell:
        """
        mkdir -p {params.outdir} logs/tobias_atacorrect
        TOBIAS ATACorrect \
            --bam {input.bam} \
            --genome {input.genome} \
            --peaks {input.peaks} \
            --blacklist {input.blacklist} \
            --outdir {params.outdir} \
            --prefix {wildcards.cond} \
            --cores {threads} > {log} 2>&1
        """


# TOBIAS ScoreBigwig: turn the bias-corrected signal into a footprint-score
# track over the consensus regions.
rule tobias_scorebigwig:
    input:
        corrected = f"{FOOTPRINT_DIR}/atacorrect/{{cond}}_corrected.bw",
        peaks     = f"{CONSENSUS_DIR}/consensus_peaks.bed"
    output:
        footprints = f"{FOOTPRINT_DIR}/{{cond}}_footprints.bw"
    threads: 16
    conda:
        "../envs/tobias.yaml"
    log:
        "logs/tobias_scorebigwig/{cond}.log"
    shell:
        """
        mkdir -p {FOOTPRINT_DIR} logs/tobias_scorebigwig
        TOBIAS ScoreBigwig \
            --signal {input.corrected} \
            --regions {input.peaks} \
            --output {output.footprints} \
            --cores {threads} > {log} 2>&1
        """


# TOBIAS BINDetect: scan the JASPAR motifs across the per-condition footprint
# tracks and report differential TF binding (ranked table + volcano + per-TF
# aggregate footprints and bound-site BEDs). `--signals` and `--cond_names`
# share CONDITION_NAMES order so they stay aligned.
rule tobias_bindetect:
    input:
        motifs  = config["jaspar_motifs"],
        signals = expand(f"{FOOTPRINT_DIR}/{{cond}}_footprints.bw", cond=CONDITION_NAMES),
        genome  = config["human_fasta"],
        peaks   = f"{CONSENSUS_DIR}/consensus_peaks.bed"
    output:
        results = f"{FOOTPRINT_DIR}/bindetect/bindetect_results.txt"
    params:
        outdir     = f"{FOOTPRINT_DIR}/bindetect",
        cond_names = " ".join(CONDITION_NAMES)
    threads: 16
    conda:
        "../envs/tobias.yaml"
    log:
        "logs/tobias_bindetect/bindetect.log"
    shell:
        """
        mkdir -p {params.outdir} logs/tobias_bindetect
        TOBIAS BINDetect \
            --motifs {input.motifs} \
            --signals {input.signals} \
            --genome {input.genome} \
            --peaks {input.peaks} \
            --cond_names {params.cond_names} \
            --outdir {params.outdir} \
            --cores {threads} > {log} 2>&1
        """
