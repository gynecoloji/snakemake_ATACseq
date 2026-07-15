# Primary ATAC-seq pipeline: FastQC/fastp → Bowtie2 alignment (human genome) →
# filtering/dedup/blacklist → MACS2 peaks → RPGC depth-normalized coverage
# → reproducible fixed-width consensus peaks + counts (Module B).
#
# Shared config, samples, directory constants and helper functions live in
# common.smk (included first by workflow/Snakefile).

# Aggregate target for the primary pipeline. Run it alone with:
#   snakemake --use-conda --cores N atacseq_all
rule atacseq_all:
    input:
        # FastQC reports
        expand(f"{FASTQC_DIR}/{{sample}}_R1_001_fastqc.html", sample=SAMPLES),
        expand(f"{FASTQC_DIR}/{{sample}}_R2_001_fastqc.html", sample=SAMPLES),
        
        # Fastp reports and trimmed reads
        expand(f"{FASTP_DIR}/{{sample}}.html", sample=SAMPLES),
        expand(f"{FASTP_DIR}/{{sample}}.json", sample=SAMPLES),

        # Raw alignment BAM (full alignment: all reads + all tags)
        expand(f"{ALIGN_DIR}/{{sample}}.bam", sample=SAMPLES),

        # Filtered BAM files
        expand(f"{FILTERED_DIR}/{{sample}}.sorted.filtered.bam", sample=SAMPLES),
        expand(f"{FILTERED_DIR}/{{sample}}.sorted.filtered.bam.bai", sample=SAMPLES),
        
        # Deduplicated BAM files
        expand(f"{DEDUP_DIR}/{{sample}}.dedup.bam", sample=SAMPLES),
        expand(f"{DEDUP_DIR}/{{sample}}.dedup.metrics.txt", sample=SAMPLES),
        
        # Blacklist filtered BAM files
        expand(f"{BLACKLIST_FILTERED_DIR}/{{sample}}.nobl.bam", sample=SAMPLES),
        expand(f"{BLACKLIST_FILTERED_DIR}/{{sample}}.nobl.bam.bai", sample=SAMPLES),
        
        # Peak calling results (with blacklist filtering)
        expand(f"{PEAKS_DIR}/{{sample}}_peaks.narrowPeak", sample=SAMPLES),
        
        # Blacklist filtering statistics
        f"{QC_DIR}/blacklist_filtering_stats.txt",

        # ── RPGC depth-normalized coverage tracks ──
        expand(f"{BIGWIG_DIR}/{{sample}}.bw", sample=SAMPLES),

        # ── Module B: consensus peaks + fragment counts ──
        f"{CONSENSUS_DIR}/consensus_peaks.bed",
        f"{CONSENSUS_DIR}/consensus_counts.txt"

# FastQC on raw reads
rule fastqc:
    input:
        r1 = "data/{sample}_R1_001.fastq.gz",
        r2 = "data/{sample}_R2_001.fastq.gz"
    output:
        r1_html = f"{FASTQC_DIR}/{{sample}}_R1_001_fastqc.html",
        r2_html = f"{FASTQC_DIR}/{{sample}}_R2_001_fastqc.html"
    params:
        outdir = FASTQC_DIR
    threads: 2
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/fastqc/{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input.r1} {input.r2} > {log} 2>&1
        """

# Fastp for read trimming and quality filtering
rule fastp:
    input:
        r1 = "data/{sample}_R1_001.fastq.gz",
        r2 = "data/{sample}_R2_001.fastq.gz"
    output:
        r1 = f"{FASTP_DIR}/{{sample}}_R1.trimmed.fastq.gz",
        r2 = f"{FASTP_DIR}/{{sample}}_R2.trimmed.fastq.gz",
        html = f"{FASTP_DIR}/{{sample}}.html",
        json = f"{FASTP_DIR}/{{sample}}.json"
    threads: 16
    conda:
        "../envs/snakemake.yaml"
    params:
        adapter_args = FASTP_ADAPTER_ARGS
    log:
        "logs/fastp/{sample}.log"
    shell:
        """
        mkdir -p {FASTP_DIR}
        fastp --in1 {input.r1} --in2 {input.r2} \
              --out1 {output.r1} --out2 {output.r2} \
              --thread {threads} \
              --html {output.html} \
              --json {output.json} \
              {params.adapter_args} \
              --trim_poly_g \
              --length_required 30 -p --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20 > {log} 2>&1
        """

# Align reads with Bowtie2 to the human index. The human FASTA must be
# chr-prefixed (UCSC) since Bowtie2 has no --add-chrname: blacklist/promoter/MACS2
# expect chr1..chrX.
rule bowtie2_align:
    input:
        r1   = f"{FASTP_DIR}/{{sample}}_R1.trimmed.fastq.gz",
        r2   = f"{FASTP_DIR}/{{sample}}_R2.trimmed.fastq.gz",
        done = f"{config['bowtie2_index']}.build.done"
    output:
        bam = f"{ALIGN_DIR}/{{sample}}.bam",
        bai = f"{ALIGN_DIR}/{{sample}}.bam.bai",
        summary = f"{ALIGN_DIR}/{{sample}}.bowtie2.log"
    params:
        index = config["bowtie2_index"]
    threads: 20
    conda:
        "../envs/snakemake.yaml"
    shell:
        """
        mkdir -p {ALIGN_DIR} {TMP_DIR}
        # Keep the FULL raw alignment (all reads incl. unmapped + multi-mappers,
        # all tags) as a coordinate-sorted, indexed BAM — nothing filtered here.
        bowtie2 -x {params.index} -1 {input.r1} -2 {input.r2} \
               -p {threads} \
               -q --phred33 -X 3000 -I 0 --no-discordant --no-mixed \
               2> {output.summary} \
            | samtools sort -@ 4 -m 2G -T {TMP_DIR}/{wildcards.sample}.rawsort -o {output.bam} -
        samtools index -@ {threads} {output.bam}
        """

# Filter (Bowtie2): properly-paired + unique HUMAN reads, drop filtering-orphaned
# mates, record per-chrom idxstats (mito-% QC), then keep only the analysis chroms
# (drop chrM/chrY/non-primary). Mirrors the CUT&RUN filter + ENCODE mito removal.
rule samtools_sort_filter_index:
    input:
        f"{ALIGN_DIR}/{{sample}}.bam"
    output:
        bam = f"{FILTERED_DIR}/{{sample}}.sorted.filtered.bam",
        bai = f"{FILTERED_DIR}/{{sample}}.sorted.filtered.bam.bai",
        flagstat = f"{FILTERED_DIR}/{{sample}}_summary.txt",
        idxstats = f"{FILTERED_DIR}/{{sample}}.idxstats.txt"
    params:
        keep_chroms    = config["keep_chroms"],
        prekeep        = f"{TMP_DIR}/{{sample}}.prekeep.bam"
    log:
        "logs/samtools/{sample}.log"
    threads: 20
    conda:
        "../envs/snakemake.yaml"
    shell:
        """
        mkdir -p {FILTERED_DIR} logs/samtools {TMP_DIR}
        # Raw flagstat of the aligned SAM (pre-filter QC)
        samtools flagstat {input} > {FILTERED_DIR}/{wildcards.sample}_raw_summary.txt 2>> {log}

        # Keep properly-paired, primary, mapped, UNIQUE reads: drop Bowtie2
        # multi-mappers (XS:i:)
        samtools view -@ {threads} -hS -f 2 -F 2316 {input} | grep -v "XS:i:" \
            > {TMP_DIR}/temp_{wildcards.sample}.unsorted.sam 2>> {log}

        # Name-sort, then drop reads orphaned by filtering (keep complete pairs only)
        samtools sort -n -O sam -o {TMP_DIR}/temp_{wildcards.sample}.sorted.sam \
            {TMP_DIR}/temp_{wildcards.sample}.unsorted.sam 2>> {log}
        python workflow/scripts/process_sam.py {TMP_DIR}/temp_{wildcards.sample}.sorted.sam \
            {TMP_DIR}/temp_{wildcards.sample}.uc.unsorted.sam 2>> {log}

        # Coordinate-sort human reads (incl chrM) to a BAM and index it
        samtools view -@ {threads} -bhS {TMP_DIR}/temp_{wildcards.sample}.uc.unsorted.sam | \
        samtools sort -@ {threads} -O bam -o {params.prekeep} 2>> {log}
        samtools index -@ {threads} {params.prekeep} 2>> {log}

        # Per-chromosome counts (mito-% QC) BEFORE dropping non-analysis chroms
        samtools idxstats {params.prekeep} > {output.idxstats} 2>> {log}

        # Keep only the analysis chromosomes (drop chrM/chrY/non-primary), index, flagstat
        samtools view -@ {threads} -b -o {output.bam} {params.prekeep} {params.keep_chroms} 2>> {log}
        samtools index -@ {threads} {output.bam} {output.bai} 2>> {log}
        samtools flagstat {output.bam} > {output.flagstat} 2>> {log}

        # Clean up temporary files
        rm -f {TMP_DIR}/temp_{wildcards.sample}.unsorted.sam \
              {TMP_DIR}/temp_{wildcards.sample}.uc.unsorted.sam \
              {TMP_DIR}/temp_{wildcards.sample}.sorted.sam \
              {params.prekeep} {params.prekeep}.bai
        """

# Remove duplicates with Picard
rule remove_duplicates:
    input:
        filtered_bam = f"{FILTERED_DIR}/{{sample}}.sorted.filtered.bam"
    output:
        dedup_bam = f"{DEDUP_DIR}/{{sample}}.dedup.bam",
        metrics = f"{DEDUP_DIR}/{{sample}}.dedup.metrics.txt"
    threads: 4
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/dedup/{sample}.log"
    shell:
        """
        mkdir -p {DEDUP_DIR}
        java -jar ref/picard.jar MarkDuplicates \
               INPUT={input.filtered_bam} \
               OUTPUT={output.dedup_bam} \
               METRICS_FILE={output.metrics} \
               REMOVE_DUPLICATES=true \
               ASSUME_SORTED=true \
               VALIDATION_STRINGENCY=LENIENT \
               TMP_DIR=tmp 2> {log}
        samtools index {output.dedup_bam}
        """

# Filter against blacklist regions
rule filter_blacklist:
    priority: 10
    input:
        bam = f"{DEDUP_DIR}/{{sample}}.dedup.bam",
        blacklist = config["blacklist"]
    output:
        filtered_bam = f"{BLACKLIST_FILTERED_DIR}/{{sample}}.nobl.bam",
        filtered_bai = f"{BLACKLIST_FILTERED_DIR}/{{sample}}.nobl.bam.bai",
        excluded_reads = f"{BLACKLIST_FILTERED_DIR}/{{sample}}.blacklisted.bam"
    params:
        temp_bedpe = f"{TMP_DIR}/{{sample}}.fragments.bedpe",
        temp_fragment_bed = f"{TMP_DIR}/{{sample}}.fragments.bed",
        temp_blacklist_fragments = f"{TMP_DIR}/{{sample}}.blacklisted.fragments.bed",
        temp_blacklist_ids = f"{TMP_DIR}/{{sample}}.blacklisted.ids.txt",
        temp_namesorted_bam = f"{TMP_DIR}/{{sample}}.namesorted.bam",
        temp_filtered_bam = f"{TMP_DIR}/{{sample}}.filtered.bam",
        temp_excluded_bam = f"{TMP_DIR}/{{sample}}.excluded.bam"
    threads: 8
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/blacklist_filter/{sample}.log"
    shell:
        """
        mkdir -p {BLACKLIST_FILTERED_DIR} {TMP_DIR}
        
        # Sort BAM by read name once and reuse for both BEDPE conversion and filtering
        samtools sort -n -@ {threads} -o {params.temp_namesorted_bam} {input.bam} 2> {log}

        # Convert name-sorted BAM to BEDPE format
        bedtools bamtobed -bedpe -i {params.temp_namesorted_bam} > {params.temp_bedpe} 2>> {log}

        # Convert BEDPE to fragment BED (one entry per fragment)
        # Extract the fragment coordinates (minimum start, maximum end)
        # and keep the read name for later filtering
        awk 'BEGIN {{OFS="\\t"}} {{if ($1==$4) print $1, ($2<$5?$2:$5), ($3>$6?$3:$6), $7, ".", ($9=="+"?"+":"-")}}' \
        {params.temp_bedpe} > {params.temp_fragment_bed} 2>> {log}

        # Find fragments that intersect with blacklisted regions
        bedtools intersect -a {params.temp_fragment_bed} -b {input.blacklist} -wa > {params.temp_blacklist_fragments} 2>> {log}

        # Extract read IDs from blacklisted fragments
        cut -f4 {params.temp_blacklist_fragments} | sort | uniq > {params.temp_blacklist_ids} 2>> {log}
        
        # Create properly paired BAMs - one with fragments that don't overlap blacklist
        
        # Filter out fragments overlapping blacklisted regions
        samtools view -@ {threads} -b -N ^{params.temp_blacklist_ids} \
            {params.temp_namesorted_bam} > {params.temp_filtered_bam} 2>> {log}
            
        # Extract fragments overlapping blacklisted regions
        samtools view -@ {threads} -b -N {params.temp_blacklist_ids} \
            {params.temp_namesorted_bam} > {params.temp_excluded_bam} 2>> {log}
            
        # Sort filtered BAM (non-blacklisted fragments) by coordinate for final output
        samtools sort -@ {threads} -o {output.filtered_bam} {params.temp_filtered_bam} 2>> {log}
        
        # Sort excluded BAM (blacklisted fragments) by coordinate for QC
        samtools sort -@ {threads} -o {output.excluded_reads} {params.temp_excluded_bam} 2>> {log}
        
        # Index the filtered BAM
        samtools index -@ {threads} {output.filtered_bam} {output.filtered_bai} 2>> {log}
        
        # Report stats (before cleanup so temp files are still available)
        echo "Blacklist filtering completed for {wildcards.sample}" >> {log}
        echo "$(wc -l < {params.temp_blacklist_fragments}) fragments overlap blacklisted regions" >> {log}
        echo "$(wc -l < {params.temp_blacklist_ids}) unique fragment IDs overlapping blacklisted regions" >> {log}
        echo "$(samtools view -c {output.excluded_reads}) total reads in excluded fragments" >> {log}
        echo "$(samtools view -c {output.filtered_bam}) total reads in filtered output" >> {log}

        # Clean up temporary files
        rm -f {params.temp_bedpe} {params.temp_fragment_bed} {params.temp_blacklist_fragments} \
            {params.temp_blacklist_ids} {params.temp_namesorted_bam} {params.temp_filtered_bam} \
            {params.temp_excluded_bam}
        """
        
# Call peaks with MACS2 (without input control) - now using blacklist filtered BAM
rule call_peaks:
    input:
        treatment = f"{BLACKLIST_FILTERED_DIR}/{{sample}}.nobl.bam"
    output:
        peaks = f"{PEAKS_DIR}/{{sample}}_peaks.narrowPeak"
    params:
        outdir = PEAKS_DIR,
        name = "{sample}"
    conda:
        "../envs/macs2.yaml"
    log:
        "logs/macs2/{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}
        macs2 callpeak \
              -t {input.treatment} \
              -f BAMPE \
              -g hs \
              --outdir {params.outdir} \
              -n {params.name} \
              --nomodel \
              -q 0.05 > {log} 2>&1
        """

# ── Build the human Bowtie2 index (once) ──
# Optionally restrict the human genome to `align_chroms` before building the index
# (reads then align only to those chroms). The human FASTA must be chr-prefixed
# (UCSC) to match the blacklist/promoter/MACS2 conventions.
rule build_genome_index:
    input:
        human = config["human_fasta"]
    output:
        done = touch(f"{config['bowtie2_index']}.build.done")
    params:
        index        = config["bowtie2_index"],
        chroms       = config["align_chroms"],
        human_subset = f"{config['bowtie2_index']}.human.subset.fa"
    threads: 8
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/build_genome_index/build.log"
    shell:
        """
        mkdir -p $(dirname {params.index}) logs/build_genome_index
        # Optionally subset the human genome to the requested chromosomes
        if [ -n "{params.chroms}" ]; then
            samtools faidx {input.human} {params.chroms} > {params.human_subset} 2>> {log}
            HUMAN={params.human_subset}
        else
            HUMAN={input.human}
        fi
        bowtie2-build --threads {threads} $HUMAN {params.index} >> {log} 2>&1
        if [ -n "{params.chroms}" ]; then rm -f {params.human_subset}; fi
        """

# ── RPGC depth-normalized bigWig ───────────────────────────────────────
rule create_bigwig:
    input:
        bam = f"{BLACKLIST_FILTERED_DIR}/{{sample}}.nobl.bam",
        bai = f"{BLACKLIST_FILTERED_DIR}/{{sample}}.nobl.bam.bai"
    output:
        bw = f"{BIGWIG_DIR}/{{sample}}.bw"
    params:
        egs       = config["effective_genome_size"],
        bin_size  = config["bin_size"],
        blacklist = config["blacklist"]
    threads: 8
    conda:
        "../envs/deeptools.yaml"
    log:
        "logs/bigwig/{sample}.log"
    shell:
        """
        mkdir -p {BIGWIG_DIR} logs/bigwig
        bamCoverage --bam {input.bam} \
            --normalizeUsing RPGC \
            --effectiveGenomeSize {params.egs} \
            --binSize {params.bin_size} \
            --numberOfProcessors {threads} \
            --extendReads \
            --blackListFileName {params.blacklist} \
            --outFileName {output.bw} > {log} 2>&1
        """

# ── Module B: relaxed MACS2 calls for IDR (2-replicate conditions only) ─
rule relaxed_peaks:
    wildcard_constraints:
        sample = _alt(IDR_SAMPLES)
    input:
        bam = f"{BLACKLIST_FILTERED_DIR}/{{sample}}.nobl.bam"
    output:
        peaks = f"{RELAXED_PEAKS_DIR}/{{sample}}_relaxed.narrowPeak"
    params:
        outdir = RELAXED_PEAKS_DIR,
        name   = "{sample}",
        genome = config["macs2_genome"],
        pvalue = config["idr_relaxed_pvalue"],
        top_n  = config["idr_top_n_peaks"]
    conda:
        "../envs/macs2.yaml"
    log:
        "logs/relaxed_peaks/{sample}.log"
    shell:
        """
        mkdir -p {params.outdir} logs/relaxed_peaks
        macs2 callpeak \
            -t {input.bam} \
            -f BAMPE -g {params.genome} \
            --outdir {params.outdir} \
            -n {params.name}_relaxedtmp \
            --nomodel -p {params.pvalue} > {log} 2>&1
        sort -k8,8gr {params.outdir}/{params.name}_relaxedtmp_peaks.narrowPeak \
            > {params.outdir}/{params.name}_relaxedtmp_sorted.narrowPeak
        head -n {params.top_n} {params.outdir}/{params.name}_relaxedtmp_sorted.narrowPeak > {output.peaks}
        rm -f {params.outdir}/{params.name}_relaxedtmp_peaks.narrowPeak \
              {params.outdir}/{params.name}_relaxedtmp_peaks.xls \
              {params.outdir}/{params.name}_relaxedtmp_summits.bed \
              {params.outdir}/{params.name}_relaxedtmp_sorted.narrowPeak
        """

# ── Module B: IDR reproducibility (exactly 2-replicate conditions) ──────
# _group_relaxed_inputs() is defined in common.smk.
rule reproducible_idr:
    wildcard_constraints:
        group = _alt(IDR_GROUPS)
    input:
        peaks = _group_relaxed_inputs
    output:
        peaks = f"{CONSENSUS_DIR}/idr/{{group}}.idr_peaks.narrowPeak"
    params:
        threshold = config["idr_threshold"],
        idr_out   = f"{CONSENSUS_DIR}/idr/{{group}}.idr.txt"
    conda:
        "../envs/idr.yaml"
    log:
        "logs/reproducible/{group}_idr.log"
    shell:
        """
        mkdir -p {CONSENSUS_DIR}/idr logs/reproducible
        idr --samples {input.peaks} \
            --input-file-type narrowPeak \
            --rank p.value \
            --idr-threshold {params.threshold} \
            --output-file {params.idr_out} > {log} 2>&1
        # Keep peaks whose global IDR (col 12, -log10) passes the threshold, and
        # re-emit standard narrowPeak columns 1-10 (score col 9 = -log10 q).
        awk -v t={params.threshold} \
            'BEGIN{{OFS="\\t"; c=-log(t)/log(10)}} $12>=c {{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' \
            {params.idr_out} > {output.peaks}
        """

# ── Module B: build the fixed-width, non-overlapping consensus set ──────
rule consensus_peaks:
    input:
        narrowpeaks = expand(f"{PEAKS_DIR}/{{sample}}_peaks.narrowPeak", sample=SAMPLES),
        idr         = expand(f"{CONSENSUS_DIR}/idr/{{group}}.idr_peaks.narrowPeak", group=IDR_GROUPS),
        blacklist   = config["blacklist"]
    output:
        bed = f"{CONSENSUS_DIR}/consensus_peaks.bed",
        saf = f"{CONSENSUS_DIR}/consensus_peaks.saf"
    params:
        groups       = GROUPS,
        group_method = GROUP_METHOD,
        peaks_dir    = PEAKS_DIR,
        idr_dir      = f"{CONSENSUS_DIR}/idr",
        min_reps     = config["consensus_min_replicates"],
        window       = config["consensus_window"],
        keep_regex   = config["keep_chroms_regex"]
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/consensus/consensus.log"
    script:
        "../scripts/consensus_peaks.py"

# ── Module B: count fragments over the consensus set (all samples) ──────
rule count_fragments_consensus:
    input:
        saf  = f"{CONSENSUS_DIR}/consensus_peaks.saf",
        bams = expand(f"{BLACKLIST_FILTERED_DIR}/{{sample}}.nobl.bam", sample=SAMPLES),
        bais = expand(f"{BLACKLIST_FILTERED_DIR}/{{sample}}.nobl.bam.bai", sample=SAMPLES)
    output:
        counts  = f"{CONSENSUS_DIR}/consensus_counts.txt",
        summary = f"{CONSENSUS_DIR}/consensus_counts.txt.summary"
    threads: 8
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/consensus_counts/featurecounts.log"
    shell:
        """
        mkdir -p {CONSENSUS_DIR} logs/consensus_counts
        featureCounts -F SAF -a {input.saf} \
            -p --countReadPairs \
            -T {threads} \
            -o {output.counts} \
            {input.bams} > {log} 2>&1
        """

# Generate blacklist filtering statistics
rule blacklist_stats:
    input:
        original_bams = expand(f"{DEDUP_DIR}/{{sample}}.dedup.bam", sample=SAMPLES),
        filtered_bams = expand(f"{BLACKLIST_FILTERED_DIR}/{{sample}}.nobl.bam", sample=SAMPLES),
        excluded_bams = expand(f"{BLACKLIST_FILTERED_DIR}/{{sample}}.blacklisted.bam", sample=SAMPLES)
    output:
        stats = f"{QC_DIR}/blacklist_filtering_stats.txt"
    params:
        samples = SAMPLES  # Pass sample names to the script
    threads: 1
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/blacklist_stats/summary.log"
    script:
        "../scripts/blacklist-stats-script.py"