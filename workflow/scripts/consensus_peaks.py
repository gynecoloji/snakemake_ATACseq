#!/usr/bin/env python3
"""Reproducible fixed-width non-overlapping consensus ATAC-seq peak set.

Corces et al. 2018 (fixed-width summit windows + score-per-million iterative
overlap removal), with replicate-count-adaptive reproducibility:
  * >=3 reps -> majority vote (summit covered by >= min_reps replicates)
  * ==2 reps -> IDR peaks (precomputed by the reproducible_idr rule)
  * ==1 rep  -> passthrough (flagged; not reproducibility-filtered)
Reproducible peaks are UNIONed across conditions, resized to fixed windows,
blacklist/chrom filtered, and reduced by SPM-ranked iterative overlap removal.
"""
import re
import bisect
from pathlib import Path
from collections import defaultdict


class Peak:
    __slots__ = ("chrom", "start", "end", "score", "summit", "sample", "spm")

    def __init__(self, chrom, start, end, score, summit, sample):
        self.chrom = chrom
        self.start = start      # original peak start (0-based)
        self.end = end          # original peak end
        self.score = score      # MACS2 -log10(q) (narrowPeak col 9)
        self.summit = summit    # absolute summit = start + col10 offset
        self.sample = sample    # sample id (or group name for IDR peaks)
        self.spm = 0.0


def load_narrowpeak(path, sample):
    """Parse a MACS2 narrowPeak file into Peak objects tagged with `sample`."""
    peaks = []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        f = line.split("\t")
        start = int(f[1])
        peaks.append(Peak(f[0], start, int(f[2]), float(f[8]), start + int(f[9]), sample))
    return peaks


def _point_covered(pos, intervals):
    """True if pos falls in any (start, end) of the sorted intervals list."""
    for s, e in intervals:
        if s <= pos < e:
            return True
        if s > pos:
            break
    return False


def majority_keep(rep_peaks, min_reps):
    """rep_peaks: list (per replicate) of Peak lists. Keep peaks whose summit is
    covered by peaks in >= min_reps replicates of the same condition."""
    rep_index = []
    for peaks in rep_peaks:
        by_chrom = defaultdict(list)
        for p in peaks:
            by_chrom[p.chrom].append((p.start, p.end))
        for c in by_chrom:
            by_chrom[c].sort()
        rep_index.append(by_chrom)

    kept = []
    for peaks in rep_peaks:
        for p in peaks:
            cover = sum(1 for idx in rep_index if _point_covered(p.summit, idx.get(p.chrom, [])))
            if cover >= min_reps:
                kept.append(p)
    return kept


def assign_spm(peaks):
    """Assign per-sample score-per-million to each peak (in place)."""
    total = defaultdict(float)
    for p in peaks:
        total[p.sample] += p.score
    for p in peaks:
        denom = total[p.sample] / 1e6
        p.spm = p.score / denom if denom > 0 else 0.0
    return peaks


def fixed_window(peak, width):
    """Return (chrom, wstart, wend): fixed-`width` window centered on the summit."""
    wstart = peak.summit - width // 2
    if wstart < 0:
        wstart = 0
    return (peak.chrom, wstart, wstart + width)


def _load_bed_intervals(path):
    by_chrom = defaultdict(list)
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.split("\t")
        by_chrom[f[0]].append((int(f[1]), int(f[2])))
    for c in by_chrom:
        by_chrom[c].sort()
    return by_chrom


def _overlaps_any(chrom, start, end, by_chrom):
    for s, e in by_chrom.get(chrom, []):
        if s < end and start < e:
            return True
        if s >= end:
            break
    return False


def iterative_overlap_removal(windows, width):
    """Greedily keep highest-SPM window; drop later windows overlapping a kept one.
    All windows are exactly `width` bp, so two overlap iff |start_a - start_b| < width."""
    windows = sorted(windows, key=lambda w: w["spm"], reverse=True)
    kept = []
    starts_by_chrom = defaultdict(list)
    for w in windows:
        starts = starts_by_chrom[w["chrom"]]
        i = bisect.bisect_left(starts, w["start"])
        overlap = (i > 0 and w["start"] - starts[i - 1] < width) or \
                  (i < len(starts) and starts[i] - w["start"] < width)
        if not overlap:
            bisect.insort(starts, w["start"])
            kept.append(w)
    return kept


def build_consensus(groups, group_method, narrowpeak_paths, idr_paths,
                    min_reps, width, keep_regex, blacklist_path):
    """Return sorted, named list of consensus windows: {chrom,start,end,name,spm}."""
    keep_re = re.compile(keep_regex)
    kept = []
    for g, members in groups.items():
        method = group_method[g]
        if method == "majority":
            reps = [load_narrowpeak(narrowpeak_paths[s], s) for s in members]
            kept.extend(majority_keep(reps, min_reps))
        elif method == "single":
            kept.extend(load_narrowpeak(narrowpeak_paths[members[0]], members[0]))
        elif method == "idr":
            kept.extend(load_narrowpeak(idr_paths[g], g))

    assign_spm(kept)

    blacklist = _load_bed_intervals(blacklist_path)
    windows = []
    for p in kept:
        if not keep_re.fullmatch(p.chrom):
            continue
        chrom, ws, we = fixed_window(p, width)
        if _overlaps_any(chrom, ws, we, blacklist):
            continue
        windows.append({"chrom": chrom, "start": ws, "end": we, "spm": p.spm})

    consensus = iterative_overlap_removal(windows, width)
    consensus.sort(key=lambda w: (w["chrom"], w["start"]))
    for i, w in enumerate(consensus, 1):
        w["name"] = f"consensus_peak_{i}"
    return consensus


def write_bed(consensus, path):
    lines = [f'{w["chrom"]}\t{w["start"]}\t{w["end"]}\t{w["name"]}\t{int(round(w["spm"]))}\t.'
             for w in consensus]
    Path(path).write_text("\n".join(lines) + ("\n" if lines else ""))


def write_saf(consensus, path):
    lines = ["GeneID\tChr\tStart\tEnd\tStrand"]
    for w in consensus:  # SAF is 1-based, inclusive
        lines.append(f'{w["name"]}\t{w["chrom"]}\t{w["start"] + 1}\t{w["end"]}\t.')
    Path(path).write_text("\n".join(lines) + "\n")


if "snakemake" in globals():  # pragma: no cover
    sm = snakemake  # noqa: F821
    _groups = dict(sm.params.groups)
    _method = dict(sm.params.group_method)
    _npaths = {s: f"{sm.params.peaks_dir}/{s}_peaks.narrowPeak"
               for m in _groups.values() for s in m}
    _ipaths = {g: f"{sm.params.idr_dir}/{g}.idr_peaks.narrowPeak"
               for g in _groups if _method[g] == "idr"}
    _consensus = build_consensus(_groups, _method, _npaths, _ipaths,
                                 int(sm.params.min_reps), int(sm.params.window),
                                 sm.params.keep_regex, str(sm.input.blacklist))
    write_bed(_consensus, str(sm.output.bed))
    write_saf(_consensus, str(sm.output.saf))
