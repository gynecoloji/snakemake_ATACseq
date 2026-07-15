#!/usr/bin/env python3
"""Generate promoter BEDs (TSS +/- N bp) under several TSS definitions.

Usage (run from the repo root):
    python workflow/scripts/build_promoter_beds.py <mode> [--windows 3000,5000]
                                             [--chroms LIST|all] [--config PATH] [--label TOKEN]
    modes: transcript | unique | mane | cage | all

TSS definitions
---------------
transcript  Every transcript's 5' end (strand-aware) from gencode v36.
            + strand -> transcript start; - strand -> transcript end.
            ~231k windows on chr1-22,X (one per transcript).  <- legacy set
unique      transcript set collapsed to distinct (chrom, TSS, strand) tuples
            (~206k) -- removes isoform double-counting, keeps every alternative
            promoter. A representative transcript supplies the annotation cols.
mane        One canonical TSS per gene: the MANE Select transcript where the
            gene has one, else the 5'-most transcript. ~one row per gene.
cage        Empirical TSS from FANTOM5 robust CAGE peaks (hg38); the peak's
            representative position (BED thickStart) is the TSS. BED6 output.

Kept chromosomes (configurable)
-------------------------------
By default the kept set is read from `keep_chroms` in config/config.yaml (so the
promoter scope tracks the analysis scope); if that is unavailable it falls back
to chr1-22,X. Override per run with --chroms (a comma-separated list, or `all`
to keep every chromosome in the reference). The output filename carries a scope
token derived from the kept set (e.g. `chr1-22X`); override it with --label.

Window: BED half-open [max(0, tss0-N), min(chromlen, tss0+N)); width up to 2N,
score column = 2N. LF line endings. GTF-derived modes emit the legacy 10-column
layout
  tx_id  chrom  start  end  score  strand  tx_id.ver  transcript_type  gene_name  gene_id
sorted by transcript_id. CAGE emits BED6 (chrom start end peak_id score strand)
sorted by genomic position.
"""
import re, sys, gzip, argparse

GTF  = "ref/gencode.v36.annotation.gtf"
FAI  = "ref/hg38.fa.fai"
CAGE = "ref/FANTOM5_CAGE_peaks_hg38.bed.gz"
DEFAULT_CHROMS = [f"chr{c}" for c in list(range(1, 23)) + ["X"]]

OUT = {
    "transcript": "ref/Promoter_UCSC_hg38_{N}bp_{label}.bed",
    "unique":     "ref/Promoter_uniqueTSS_hg38_{N}bp_{label}.bed",
    "mane":       "ref/Promoter_MANEcanonical_hg38_{N}bp_{label}.bed",
    "cage":       "ref/Promoter_FANTOM5CAGE_hg38_{N}bp_{label}.bed",
}


def read_config_chroms(path):
    """Pull keep_chroms (a YAML flow list) from config.yaml; None if not found."""
    try:
        txt = open(path).read()
    except OSError:
        return None
    m = re.search(r'(?m)^\s*keep_chroms:\s*\[([^\]]*)\]', txt)
    if not m:
        return None
    chroms = [c.strip().strip('"\'') for c in m.group(1).split(",") if c.strip()]
    return chroms or None


def chrom_label(chroms):
    """Compact scope token, e.g. {chr1..chr22,chrX} -> 'chr1-22X'."""
    bare = lambda c: c[3:] if c.startswith("chr") else c
    nums = sorted(int(bare(c)) for c in chroms if bare(c).isdigit())
    letters = sorted(bare(c) for c in chroms if not bare(c).isdigit())
    parts, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        parts.append(str(nums[i]) if i == j else f"{nums[i]}-{nums[j]}")
        i = j + 1
    return "chr" + "".join(parts) + "".join(letters)


def chrom_sizes(keep):
    """{chrom: length} from the .fai, restricted to `keep` (None = all chroms)."""
    size = {}
    with open(FAI) as f:
        for ln in f:
            p = ln.split("\t")
            if keep is None or p[0] in keep:
                size[p[0]] = int(p[1])
    return size


def parse_gtf_transcripts(size):
    rx_tid = re.compile(r'transcript_id "([^"]+)"')
    rx_tt  = re.compile(r'transcript_type "([^"]+)"')
    rx_gn  = re.compile(r'gene_name "([^"]+)"')
    rx_gid = re.compile(r'gene_id "([^"]+)"')
    recs = []
    with open(GTF) as f:
        for line in f:
            if "\ttranscript\t" not in line:
                continue
            c = line.rstrip("\n").split("\t")
            if c[2] != "transcript" or c[0] not in size:
                continue
            a = c[8]
            tid_v = rx_tid.search(a).group(1); tid = tid_v.split(".")[0]
            m = rx_tt.search(a);  tt = m.group(1) if m else ""
            m = rx_gn.search(a);  gn = m.group(1) if m else ""
            m = rx_gid.search(a); gid = m.group(1).split(".")[0] if m else ""
            strand = c[6]
            tss0 = (int(c[3]) - 1) if strand == "+" else (int(c[4]) - 1)
            recs.append({"chrom": c[0], "tss0": tss0, "strand": strand,
                         "tid_v": tid_v, "tid": tid, "tt": tt, "gn": gn, "gid": gid,
                         "mane": 'tag "MANE_Select"' in a})
    return recs


def select(mode, recs):
    # Where several transcripts collapse to one row (unique TSS, or a gene's
    # canonical TSS), the representative is chosen by a stable transcript_id
    # tie-break so the output is reproducible regardless of GTF line order.
    if mode == "transcript":
        return recs
    if mode == "unique":
        best = {}
        for r in recs:
            k = (r["chrom"], r["tss0"], r["strand"])
            if k not in best or r["tid"] < best[k]["tid"]:
                best[k] = r
        return list(best.values())
    if mode == "mane":
        by_gene = {}
        for r in recs:
            by_gene.setdefault(r["gid"], []).append(r)
        out = []
        for gid, rs in by_gene.items():
            mane = [r for r in rs if r["mane"]]
            if mane:                                            # MANE Select (stable if >1)
                out.append(min(mane, key=lambda r: r["tid"]))
            else:                                               # 5'-most transcript, tie -> lowest tid
                plus = rs[0]["strand"] == "+"
                key = (lambda r: (r["tss0"], r["tid"])) if plus else (lambda r: (-r["tss0"], r["tid"]))
                out.append(min(rs, key=key))
        return out
    raise ValueError(mode)


def write_gtf_mode(mode, recs, windows, size, label):
    rows = select(mode, recs)
    rows.sort(key=lambda r: r["tid"])
    for N in windows:
        path = OUT[mode].format(N=N, label=label)
        with open(path, "w", newline="\n") as o:
            for r in rows:
                s = max(0, r["tss0"] - N); e = min(size[r["chrom"]], r["tss0"] + N)
                o.write("\t".join(map(str, [r["tid"], r["chrom"], s, e, 2 * N,
                        r["strand"], r["tid_v"], r["tt"], r["gn"], r["gid"]])) + "\n")
        print(f"  wrote {path}: {len(rows)} rows")


def write_cage(windows, size, label):
    peaks = []
    with gzip.open(CAGE, "rt") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 8 or c[0] not in size:
                continue
            peaks.append((c[0], int(c[6]), c[3], c[4], c[5]))   # chrom, thickStart(TSS0), name, score, strand
    peaks.sort(key=lambda p: (p[0], p[1]))
    for N in windows:
        path = OUT["cage"].format(N=N, label=label)
        with open(path, "w", newline="\n") as o:
            for chrom, tss0, name, score, strand in peaks:
                s = max(0, tss0 - N); e = min(size[chrom], tss0 + N)
                o.write("\t".join([chrom, str(s), str(e), name, str(score), strand]) + "\n")
        print(f"  wrote {path}: {len(peaks)} rows")


def main(argv):
    ap = argparse.ArgumentParser(description="Generate promoter BEDs under several TSS definitions.")
    ap.add_argument("mode", choices=["transcript", "unique", "mane", "cage", "all"])
    ap.add_argument("--windows", default="3000,5000", help="comma-separated half-widths (bp)")
    ap.add_argument("--chroms", default=None,
                    help="chromosomes to keep: comma-separated list or 'all'. "
                         "Default: keep_chroms from --config, else chr1-22,X.")
    ap.add_argument("--config", default="config/config.yaml",
                    help="config YAML providing the default keep_chroms (default: config/config.yaml)")
    ap.add_argument("--label", default=None,
                    help="filename scope token (default: derived from the kept set, e.g. chr1-22X)")
    a = ap.parse_args(argv)

    windows = [int(x) for x in a.windows.split(",")]
    if a.chroms:
        keep = None if a.chroms.strip().lower() == "all" \
               else [c.strip() for c in a.chroms.split(",") if c.strip()]
    else:
        keep = read_config_chroms(a.config) or DEFAULT_CHROMS
    keep_set = set(keep) if keep is not None else None
    size = chrom_sizes(keep_set)
    label = a.label or ("allchr" if keep_set is None else chrom_label(size.keys()))
    print(f"keep chroms: {label} ({len(size)} chromosomes); windows: {windows}")

    modes = ["transcript", "unique", "mane", "cage"] if a.mode == "all" else [a.mode]
    recs = parse_gtf_transcripts(size) if any(m != "cage" for m in modes) else None
    for m in modes:
        print(f"[{m}] TSS definition:")
        if m == "cage":
            write_cage(windows, size, label)
        else:
            write_gtf_mode(m, recs, windows, size, label)


if __name__ == "__main__":
    main(sys.argv[1:])
