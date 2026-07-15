#!/usr/bin/env python3
"""Numeric TSS enrichment score per sample from a deepTools plotProfile data table.

TSS enrichment = signal at the TSS (center bin) / mean signal in the flanking
background (outer edge bins). ~1 means no enrichment; higher is better (ENCODE
uses this as a key ATAC signal-to-noise metric).
"""
from pathlib import Path


def enrichment(values, edge_frac=0.1):
    """Center-bin value / mean of the outer `edge_frac` bins on each side."""
    n = len(values)
    if n == 0:
        return 0.0
    center = values[n // 2]
    k = max(1, int(n * edge_frac))
    bg = values[:k] + values[-k:]
    m = sum(bg) / len(bg)
    return center / m if m > 0 else 0.0


def parse_profile(path):
    """Parse `plotProfile --outFileNameData` into {sample: [profile floats]}.

    Header row (bin positions) is skipped; each data row starts with label
    field(s) followed by the numeric per-bin mean profile — we take the row's
    first field as the sample and every field that parses as a float as signal.
    """
    out = {}
    lines = [ln for ln in Path(path).read_text().splitlines() if ln.strip()]
    for line in lines[1:]:
        f = line.split("\t")
        vals = []
        for x in f[1:]:
            try:
                vals.append(float(x))
            except ValueError:
                pass
        if vals:
            out[f[0]] = vals
    return out


def build(profile_path):
    prof = parse_profile(profile_path)
    return sorted((s, round(enrichment(v), 3)) for s, v in prof.items())


def write_tsv(rows, path):
    lines = ["sample\ttss_enrichment"]
    lines += [f"{s}\t{score}" for s, score in rows]
    Path(path).write_text("\n".join(lines) + "\n")


def write_mqc(rows, path):
    lines = [
        "# id: tss_enrichment",
        "# section_name: 'TSS enrichment'",
        "# description: 'TSS enrichment score (center/background). Higher is better; <5 is poor for ATAC.'",
        "# plot_type: 'bargraph'",
        "# pconfig:",
        "#    id: 'tss_enrichment_plot'",
        "#    title: 'TSS enrichment'",
        "sample\tTSS enrichment",
    ]
    lines += [f"{s}\t{score}" for s, score in rows]
    Path(path).write_text("\n".join(lines) + "\n")


if "snakemake" in globals():  # pragma: no cover
    _rows = build(str(snakemake.input.profile))  # noqa: F821
    write_tsv(_rows, str(snakemake.output.tsv))  # noqa: F821
    write_mqc(_rows, str(snakemake.output.mqc))  # noqa: F821
