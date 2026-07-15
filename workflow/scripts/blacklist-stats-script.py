#!/usr/bin/env python3
"""
Generate statistics about reads filtered with blacklist regions.
"""
import pandas as pd
import os
import sys
from snakemake.shell import shell

# Function to count reads in a BAM file
def count_reads(bam_file):
    cmd = f"samtools view -c {bam_file}"
    return int(shell(cmd, read=True).strip())

# Access input and output from Snakemake
original_bam_files = snakemake.input.original_bams
filtered_bam_files = snakemake.input.filtered_bams
excluded_bam_files = snakemake.input.excluded_bams
stats_file = snakemake.output.stats

# Get sample names
samples = snakemake.params.get("samples", [os.path.basename(f).replace(".dedup.bam", "") for f in original_bam_files])

stats = []

# Create stats directory if it doesn't exist
os.makedirs(os.path.dirname(stats_file), exist_ok=True)

# Process each sample
for i, sample in enumerate(samples):
    # Count reads
    original_count = count_reads(original_bam_files[i])
    filtered_count = count_reads(filtered_bam_files[i])
    excluded_count = count_reads(excluded_bam_files[i])
    
    # Calculate percentage
    percent_excluded = (excluded_count / original_count * 100) if original_count > 0 else 0
    
    stats.append({
        "Sample": sample,
        "Original_Reads": original_count,
        "Filtered_Reads": filtered_count,
        "Blacklisted_Reads": excluded_count,
        "Percent_Excluded": f"{percent_excluded:.2f}%"
    })

# Create DataFrame and write to file
df = pd.DataFrame(stats)
with open(stats_file, "w") as f:
    f.write("# Blacklist Filtering Statistics\n")
    f.write("# Date: " + shell("date", read=True).strip() + "\n\n")
    f.write(df.to_string(index=False))
    f.write("\n\nTotal reads before filtering: " + str(df["Original_Reads"].sum()))
    f.write("\nTotal reads after filtering: " + str(df["Filtered_Reads"].sum()))
    f.write("\nTotal blacklisted reads: " + str(df["Blacklisted_Reads"].sum()))
    f.write("\nAverage percentage excluded: " + 
            f"{df['Blacklisted_Reads'].sum() / df['Original_Reads'].sum() * 100:.2f}%\n")
