#!/usr/bin/env python3
"""Build a portable *constitutive* CTCF annotation for hg38/GRCh38.

Why
---
The ENCODE SCREEN "CTCF-only" cCRE set (`ref/GRCh38-cCREs.CTCF-only.bed`) marks
*candidate* CTCF elements, but most of them are cell-type **variable** — only a
minority are bound across many cell types. Normalizing ATAC signal to a variable
anchor set defeats the purpose. This script grades every candidate cCRE by how
many independent ENCODE cell types have a CTCF ChIP-seq peak there, and keeps the
ones bound in >= `--min-frac` of them.

The output is a plain genome-coordinate BED, derived only from ENCODE reference
data -- it is **not** tied to any particular ATAC dataset's peak set, so the same
file can be reused for any human sample.

Method (mirrors Wang et al. 2012-style occupancy grading)
--------------------------------------------------------
1. Discover CTCF TF ChIP-seq "IDR thresholded peaks" on GRCh38 via the ENCODE
   REST API (`preferred_default=true` -> one canonical file per experiment).
   NOTE the API gotcha: `field=a,b,c` returns nothing; `field` must be REPEATED.
2. Deduplicate to one file per biosample, then take an evenly-spaced spread of
   `--n-celltypes` cell types (deterministic, so runs are reproducible).
3. Download those narrowPeak.gz files (cached; re-runs skip existing files).
4. For each candidate cCRE, count how many cell types have an overlapping CTCF
   peak, using a per-chromosome vectorized interval sweep (cumulative-max trick).
5. Threshold on occupancy fraction and write a BED6 plus the full occupancy TSV
   (so other thresholds can be re-derived later WITHOUT re-downloading).

Usage
-----
    python workflow/scripts/build_constitutive_ctcf.py \\
        --ccre ref/GRCh38-cCREs.CTCF-only.bed \\
        --out  ref/constitutive_ctcf_hg38.bed \\
        [--occupancy ref/ctcf_occupancy_hg38.tsv] \\
        [--n-celltypes 60] [--min-frac 0.90] [--cache <dir>] [--jobs 8]
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

ENCODE = "https://www.encodeproject.org"
SEARCH = (
    ENCODE + "/search/?type=File"
    "&assembly=GRCh38"
    "&file_type=bed+narrowPeak"
    "&assay_title=TF+ChIP-seq"
    "&target.label=CTCF"
    "&preferred_default=true"
    "&status=released"
    # `field` MUST be repeated -- a comma-separated list silently returns nothing.
    "&field=href&field=output_type"
    "&field=biosample_ontology.term_name"
    "&field=biosample_ontology.classification"
    "&limit=all&format=json"
)


def fetch_json(url: str, timeout: int = 180) -> dict:
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)


def discover_files() -> pd.DataFrame:
    """One canonical IDR-thresholded CTCF peak file per biosample."""
    graph = fetch_json(SEARCH).get("@graph", [])
    rows = []
    for g in graph:
        if g.get("output_type") != "IDR thresholded peaks":
            continue
        bo = g.get("biosample_ontology") or {}
        href = g.get("href")
        if not href or not bo.get("term_name"):
            continue
        rows.append({
            "href": ENCODE + href,
            "accession": os.path.basename(href).split(".")[0],
            "biosample": bo["term_name"],
            "classification": bo.get("classification", ""),
        })
    df = pd.DataFrame(rows)
    if df.empty:
        sys.exit("ENCODE query returned no IDR-thresholded CTCF files")
    # one file per biosample (first hit is deterministic after sorting)
    df = df.sort_values(["biosample", "accession"]).drop_duplicates("biosample")
    return df.reset_index(drop=True)


def select_spread(df: pd.DataFrame, n: int) -> pd.DataFrame:
    """Evenly-spaced subset across the available cell types (deterministic)."""
    if n >= len(df):
        return df
    idx = np.linspace(0, len(df) - 1, n).round().astype(int)
    return df.iloc[np.unique(idx)].reset_index(drop=True)


def download(row, cache: str) -> str | None:
    dest = os.path.join(cache, f"{row.accession}.bed.gz")
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    tmp = dest + ".part"
    for attempt in range(3):
        try:
            urllib.request.urlretrieve(row.href, tmp)
            os.replace(tmp, dest)
            return dest
        except Exception as e:  # noqa: BLE001 - network flakiness is expected
            print(f"  retry {attempt+1} {row.accession}: {e}", file=sys.stderr)
    print(f"  FAILED {row.accession}", file=sys.stderr)
    return None


def read_peaks(path: str) -> pd.DataFrame:
    with gzip.open(path, "rt") as fh:
        return pd.read_csv(fh, sep="\t", header=None, usecols=[0, 1, 2],
                           names=["chrom", "start", "end"])


def build_union(paths: list[str], min_peaks: int) -> tuple[pd.DataFrame, list[str]]:
    """Union of CTCF peaks across cell types (merged), + the files actually used.

    Candidates are derived from the ChIP data itself rather than from a
    pre-classified cCRE set: SCREEN's `CA-CTCF` class deliberately excludes CTCF
    sites at promoters/enhancers, so it misses most real CTCF binding.
    Files with implausibly few peaks (failed/degenerate experiments) are dropped,
    since they can only ever depress the occupancy fraction.
    """
    frames, used = [], []
    for p in paths:
        df = read_peaks(p)
        if len(df) < min_peaks:
            print(f"      skip {os.path.basename(p)}: only {len(df)} peaks", flush=True)
            continue
        frames.append(df)
        used.append(p)
    allp = pd.concat(frames, ignore_index=True)

    blocks = []
    for chrom, sub in allp.groupby("chrom", sort=False):
        s = sub["start"].to_numpy(np.int64)
        e = sub["end"].to_numpy(np.int64)
        order = np.argsort(s, kind="stable")
        s, e = s[order], e[order]
        run_max = np.maximum.accumulate(e)
        new = np.empty(len(s), dtype=bool)
        new[0] = True
        new[1:] = s[1:] > run_max[:-1]            # gap from everything so far -> new block
        heads = np.flatnonzero(new)
        blocks.append(pd.DataFrame({
            "chrom": chrom,
            "start": s[heads],
            "end": np.maximum.reduceat(e, heads),
        }))
    union = pd.concat(blocks, ignore_index=True).sort_values(
        ["chrom", "start"]).reset_index(drop=True)
    return union, used


def peaks_by_chrom(path: str) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """narrowPeak.gz -> {chrom: (starts_sorted, cummax_ends)}."""
    with gzip.open(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t", header=None, usecols=[0, 1, 2],
                         names=["chrom", "start", "end"])
    out = {}
    for chrom, sub in df.groupby("chrom", sort=False):
        s = sub["start"].to_numpy(np.int64)
        e = sub["end"].to_numpy(np.int64)
        order = np.argsort(s, kind="stable")
        s = s[order]
        out[chrom] = (s, np.maximum.accumulate(e[order]))
    return out


def overlaps(cand_start, cand_end, starts, cummax_end) -> np.ndarray:
    """Half-open overlap: peak.start < region.end AND peak.end > region.start."""
    k = np.searchsorted(starts, cand_end, side="left")   # peaks starting before region end
    hit = np.zeros(len(cand_start), dtype=bool)
    nz = k >= 1
    hit[nz] = cummax_end[k[nz] - 1] > cand_start[nz]     # ...any ending after region start
    return hit


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ccre", default="ref/GRCh38-cCREs.CTCF-only.bed",
                    help="candidate CTCF cCRE BED (genome-wide)")
    ap.add_argument("--out", default="ref/constitutive_ctcf_hg38.bed")
    ap.add_argument("--occupancy", default="ref/ctcf_occupancy_hg38.tsv",
                    help="full per-candidate occupancy table (re-threshold offline)")
    ap.add_argument("--n-celltypes", type=int, default=60,
                    help="how many distinct cell types to grade against (default 60; "
                         "well past saturation for 'constitutive')")
    ap.add_argument("--min-frac", type=float, default=0.90,
                    help="keep cCREs bound in >= this fraction of cell types")
    ap.add_argument("--cache", default="ref/encode_ctcf",
                    help="directory for downloaded ENCODE peak files")
    ap.add_argument("--min-peaks", type=int, default=5000,
                    help="drop ENCODE files with fewer peaks than this (failed/degenerate "
                         "experiments would only depress the occupancy fraction)")
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args()

    os.makedirs(a.cache, exist_ok=True)

    print("[1/4] discovering ENCODE CTCF files ...", flush=True)
    files = discover_files()
    print(f"      {len(files)} distinct cell types available")
    files = select_spread(files, a.n_celltypes)
    print(f"      grading against {len(files)} cell types")

    print("[2/4] downloading peak files ...", flush=True)
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        paths = list(ex.map(lambda r: download(r, a.cache),
                            [r for r in files.itertuples()]))
    paths = [p for p in paths if p]
    if not paths:
        sys.exit("no ENCODE peak files could be downloaded")
    print(f"      {len(paths)} files ready")

    print("[3/4] building union CTCF candidate set ...", flush=True)
    cand, paths = build_union(paths, a.min_peaks)
    cand["name"] = "CTCF_" + cand.index.astype(str)
    w = cand["end"] - cand["start"]
    print(f"      {len(cand):,} union regions from {len(paths)} usable cell types")
    print(f"      width: median {int(w.median())} bp | p99 {int(w.quantile(0.99))} bp "
          f"| max {int(w.max())} bp")

    cs = cand["start"].to_numpy(np.int64)
    ce = cand["end"].to_numpy(np.int64)
    by_chrom = {c: sub.index.to_numpy() for c, sub in cand.groupby("chrom", sort=False)}

    counts = np.zeros(len(cand), dtype=np.int32)
    for i, p in enumerate(paths, 1):
        pk = peaks_by_chrom(p)
        hit = np.zeros(len(cand), dtype=bool)
        for chrom, idx in by_chrom.items():
            got = pk.get(chrom)
            if got is None:
                continue
            starts, cummax_end = got
            hit[idx] = overlaps(cs[idx], ce[idx], starts, cummax_end)
        counts += hit
        if i % 10 == 0 or i == len(paths):
            print(f"      {i}/{len(paths)} files", flush=True)

    frac = counts / len(paths)
    cand["ctcf_occupancy"] = counts
    cand["ctcf_frac"] = frac.round(4)

    print("[4/4] writing outputs ...", flush=True)
    cand.to_csv(a.occupancy, sep="\t", index=False)

    keep = cand[frac >= a.min_frac].copy()
    keep["score"] = np.minimum(1000, (keep["ctcf_frac"] * 1000).round().astype(int))
    keep["strand"] = "."
    keep[["chrom", "start", "end", "name", "score", "strand"]].to_csv(
        a.out, sep="\t", header=False, index=False)

    print(f"\nunion CTCF regions   : {len(cand):,}")
    print(f"cell types graded    : {len(paths)}")
    for t in (0.95, 0.90, 0.80, 0.50):
        print(f"  >= {t:.0%} occupancy : {int((frac >= t).sum()):,}")
    print(f"\nconstitutive (>= {a.min_frac:.0%}) -> {a.out}  ({len(keep):,} regions)")
    print(f"full occupancy table  -> {a.occupancy}")


if __name__ == "__main__":
    main()
