#!/usr/bin/env python3
"""Build a single self-contained interactive HTML QC report for the ATAC-seq
pipeline (everything except FastQC, which stays in MultiQC).

Stdlib only. Parse pipeline outputs into a DATA dict (interactive chart series,
not embedded images), emit HEAD+<script>DATA</script>+JS.
"""
import argparse, csv, glob, json, os, re

# Colorblind-safe (Paul Tol) qualitative palette; assigned to samples by order so
# a given sample keeps one color across every figure.
PALETTE = ["#4477aa", "#66ccee", "#228833", "#ccbb44", "#ee6677", "#aa3377",
           "#000000", "#bbbbbb", "#ee8866", "#77aadd", "#44bb99", "#aaaa00"]


def sample_colors(samples):
    """Stable sample -> hex color by position (cycles if > len(PALETTE))."""
    return {s: PALETTE[i % len(PALETTE)] for i, s in enumerate(samples)}


def parse_bowtie2_log(text):
    """Bowtie2 stderr log -> overall alignment rate + concordant-unique % + raw pairs.
    For paired-end input the leading 'N reads' line counts read PAIRS."""
    overall = re.search(r"([\d.]+)% overall alignment rate", text)
    conc = re.search(r"\(([\d.]+)%\) aligned concordantly exactly 1 time", text)
    total = re.search(r"(\d+)\s+reads; of these", text)
    return {"overall_rate": float(overall.group(1)) if overall else None,
            "concordant_uniq_pct": float(conc.group(1)) if conc else None,
            "total_pairs": int(total.group(1)) if total else None}


def mito_pct_from_idxstats(text):
    """samtools idxstats text -> chrM mapped reads as % of all mapped reads."""
    total = mito = 0
    for ln in text.splitlines():
        parts = ln.split("\t")
        if len(parts) < 3:
            continue
        chrom, mapped = parts[0], int(parts[2])
        total += mapped
        if chrom in ("chrM", "MT", "chrMT"):
            mito += mapped
    return (100.0 * mito / total) if total else 0.0


# Reads the first LIBRARY metrics block (one library per sample in this pipeline).
def parse_picard_dup(text):
    """Picard MarkDuplicates metrics text -> PERCENT_DUPLICATION as percent (0-100)."""
    lines = text.splitlines()
    for i, ln in enumerate(lines):
        if ln.startswith("LIBRARY\t") and "PERCENT_DUPLICATION" in ln:
            cols = ln.split("\t")
            idx = cols.index("PERCENT_DUPLICATION")
            for vals in lines[i + 1:]:
                if vals.strip():
                    return float(vals.split("\t")[idx]) * 100.0
    return None


def parse_complexity(text):
    """*_complexity.txt -> {NRF, PBC1, PBC2}."""
    out = {}
    for ln in text.splitlines():
        m = re.match(r"(NRF|PBC1|PBC2)\b[^\t]*\t\s*([\d.]+)", ln)
        if m:
            out[m.group(1)] = float(m.group(2))
    return out


def parse_blacklist_stats(text):
    """qc/blacklist_filtering_stats.txt (whitespace-aligned text dump) -> list of
    per-sample dicts with sample/original_reads/filtered_reads/blacklisted_reads/pct_excluded."""
    rows = []
    in_table = False
    for ln in text.splitlines():
        if not in_table:
            if "Sample" in ln and "Percent_Excluded" in ln:
                in_table = True
            continue
        if not ln.strip() or ln.strip().startswith("Total"):
            break
        fields = ln.split()
        if len(fields) == 5:
            rows.append({
                "sample": fields[0],
                "original_reads": fields[1],
                "filtered_reads": fields[2],
                "blacklisted_reads": fields[3],
                "pct_excluded": fields[4],
            })
    return rows


# --- deepTools data tables -> chart series (color-agnostic; colors added in build_data) ---
TSS_UPSTREAM = 2000    # matches computeMatrix -b/-a in qc.smk (deeptools_tss)
TSS_DOWNSTREAM = 2000


def parse_tss_profile(text):
    """Profile_TSS.data.tab -> per-sample TSS coverage curves.
    Line 0 = 'bin labels', line 1 = 'bins ...', lines 2+ = '<sample>\\tgenes\\t<v1>..<vN>'."""
    lines = [ln for ln in text.splitlines() if ln.strip()]
    series = []
    for ln in lines[2:]:
        f = ln.split("\t")
        vals = [float(x) for x in f[2:] if x != ""]
        n = len(vals)
        span = TSS_UPSTREAM + TSS_DOWNSTREAM
        pts = [[-TSS_UPSTREAM + (i + 0.5) * span / n, v] for i, v in enumerate(vals)]
        series.append({"sample": f[0], "points": pts})
    return {"x_label": "Distance from TSS (bp)", "y_label": "Mean coverage", "series": series}


def parse_gc_bias(text):
    """<sample>.gc_content.txt (computeGCBias -freq) -> [[GC%, bias_ratio], ...].
    Column 3 is the plot-ready observed/expected ratio (1.0 = unbiased)."""
    rows = [ln.split() for ln in text.splitlines() if ln.strip()]
    n = len(rows)
    return [[(i / (n - 1) * 100.0 if n > 1 else 0.0), float(r[2])] for i, r in enumerate(rows)]


def parse_pca(text):
    """deeptools_PCA.tab -> {pc1_var, pc2_var, points:[{sample,x,y}]}.
    Header lists Component + sample cols; each data row = comp, <coord per sample>, eigenvalue."""
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    header = lines[0].split("\t")
    labels = header[1:]
    if labels and labels[-1].strip().lower() == "eigenvalue":
        labels = labels[:-1]                       # trailing 'Eigenvalue' column label, not a sample
    samples = [re.sub(r"\.nobl\.bam$", "", s) for s in labels]
    ns = len(samples)
    coords, eigs = {}, []
    for ln in lines[1:]:
        f = ln.split("\t")
        comp = f[0]
        coords[comp] = [float(x) for x in f[1:1 + ns]]
        if len(f) > 1 + ns:
            eigs.append(float(f[1 + ns]))
    total = sum(eigs) or 1.0
    pc1 = coords.get("1", [0.0] * ns)
    pc2 = coords.get("2", [0.0] * ns)
    points = [{"sample": samples[i], "x": pc1[i], "y": pc2[i]} for i in range(ns)]
    pc1_var = round(100 * eigs[0] / total, 1) if eigs else None
    pc2_var = round(100 * eigs[1] / total, 1) if len(eigs) > 1 else None
    return {"pc1_var": pc1_var, "pc2_var": pc2_var, "points": points}


def parse_fingerprint(text, n_points=500):
    """ATACseq_fingerprint.tab (--outRawCounts) -> per-sample cumulative fingerprint curves.
    Line 0 = comment, line 1 = quoted bam paths, rows = per-sample bin counts. For each
    sample: sort counts asc, cumulative-sum, normalize to 1, downsample to n_points."""
    lines = text.splitlines()
    hdr = [h.strip().strip("'\"") for h in lines[1].split("\t")]
    samples = [re.sub(r"\.nobl\.bam$", "", os.path.basename(h)) for h in hdr]
    cols = [[] for _ in samples]
    for ln in lines[2:]:
        if not ln.strip():
            continue
        f = ln.split("\t")
        for j in range(len(samples)):
            cols[j].append(float(f[j]))
    series = []
    for j, s in enumerate(samples):
        vals = sorted(cols[j])
        n = len(vals)
        total = sum(vals) or 1.0
        cum, run, pts = 0.0, 0.0, []
        # pick <=n_points evenly spaced ranks; ensure last point always included
        step = max(1, (n - 1) // (n_points - 1) + 1) if n_points > 1 else n
        run = 0.0
        for k in range(n):
            run += vals[k]
            if k % step == 0 or k == n - 1:
                pts.append([k / (n - 1) if n > 1 else 0.0, run / total])
        series.append({"sample": s, "points": pts})
    return {"x_label": "Rank (fraction of bins)",
            "y_label": "Cumulative fraction of reads", "series": series}


def parse_fragment_lengths(text):
    """fragment_lengths.txt (bamPEFragmentSize --outRawFragmentLengths) -> per-sample
    (length, count) distribution. Column indices are read from the header row."""
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    head = [h.strip().lower() for h in lines[0].split("\t")]
    def col(name):
        return next(i for i, h in enumerate(head) if name in h)
    ci_size, ci_occ, ci_samp = col("size"), col("occurrence"), col("sample")
    by = {}
    order = []
    for ln in lines[1:]:
        f = ln.split("\t")
        s = re.sub(r"\.nobl\.bam$", "", os.path.basename(f[ci_samp].strip().strip("'\"")))
        if s not in by:
            by[s] = []
            order.append(s)
        by[s].append([float(f[ci_size]), float(f[ci_occ])])
    series = [{"sample": s, "points": sorted(by[s])} for s in order]
    return {"mode": "dist", "x_label": "Fragment length (bp)",
            "y_label": "Count", "series": series}


def parse_fragment_summary(text):
    """fragmentsize.txt (bamPEFragmentSize --table) -> per-sample min/q1/median/q3/max
    summary (fallback when the raw-lengths file is absent)."""
    lines = [ln for ln in text.splitlines() if ln.strip()]
    head = lines[0].split("\t")
    def col(sub):
        return next(i for i, h in enumerate(head) if sub in h)
    ci = {"min": col("Min."), "q1": col("1st"), "median": col("Median"),
          "q3": col("3rd"), "max": col("Max")}
    rows = []
    for ln in lines[1:]:
        f = ln.split("\t")
        s = re.sub(r"\.nobl\.bam$", "", os.path.basename(f[0].strip()))
        rows.append({"sample": s, **{k: float(f[i]) for k, i in ci.items()}})
    return {"mode": "summary", "rows": rows}


def nucleosome_stats(points, shoulder=12):
    """Fragment-length histogram [[len,count],...] -> ENCODE nucleosome QC.
    NFR < 147 bp; mononucleosome 147-294 bp (147..2x147); di 294-441 bp. Good ATAC
    shows a sub-nucleosomal (NFR) peak AND a mononucleosome peak. Returns the three
    band fractions, the detected mono-peak position, and pass/fail flags for each peak.

    The NFR peak is the global mode of the distribution (must be sub-nucleosomal,
    < 147 bp). The mononucleosome peak is detected as a turning point (local maximum)
    in the 150-260 bp band: a genuine bump/shoulder that rises above the decaying tail
    on both sides by `shoulder` bp. This is robust to the monotonic decay that hides
    the mononucleosome as a shoulder rather than a clean isolated peak."""
    cnt, total = {}, 0.0
    for L, c in points:
        Li, cv = int(round(float(L))), float(c)
        cnt[Li] = cnt.get(Li, 0.0) + cv
        total += cv
    if total <= 0:
        return {"nfr_pct": None, "mono_pct": None, "di_pct": None, "mono_peak_bp": None,
                "nfr_peak": "na", "mono_nucleosome": "na"}
    fr = lambda lo, hi: round(100.0 * sum(v for L, v in cnt.items() if lo <= L < hi) / total, 1)
    N, half = 520, 5                                   # 11-bp moving average over 0..520 bp
    arr = [cnt.get(L, 0.0) for L in range(N + 1)]
    sm = [sum(arr[max(0, i-half):i+half+1]) / (min(N, i+half) - max(0, i-half) + 1) for i in range(N + 1)]
    mode = max(range(N + 1), key=lambda i: sm[i])      # global peak of the distribution (NFR)
    mono, best = None, 0.0                              # mononucleosome turning point
    for i in range(150, 261):
        prom = min(sm[i] - sm[i - shoulder], sm[i] - sm[i + shoulder])
        if prom > 0 and prom >= 0.02 * sm[i] and prom > best:
            best, mono = prom, i
    return {"nfr_pct": fr(0, 147), "mono_pct": fr(147, 294), "di_pct": fr(294, 441),
            "mono_peak_bp": mono,
            "nfr_peak": "pass" if mode < 147 else "warn",
            "mono_nucleosome": "pass" if mono is not None else "fail"}


def parse_cor_matrix(text):
    """correlation_matrix.tab (plotCorrelation --outFileCorMatrix) -> {labels, matrix}.
    Header = leading tab + quoted sample labels; each row = quoted label + values."""
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    def clean(x):
        return re.sub(r"\.nobl\.bam$", "", x.strip().strip("'\""))
    labels = [clean(x) for x in lines[0].split("\t")[1:]]
    matrix = [[float(v) for v in ln.split("\t")[1:]] for ln in lines[1:]]
    return {"labels": labels, "matrix": matrix}


def parse_tss_heatmap(text):
    """tss_heatmap_downsampled.json (from workflow/scripts/downsample_tss_matrix.py) -> dict, passthrough."""
    return json.loads(text)


def read_tsv(path):
    """Read a tab-separated file with a header row into a list of dicts."""
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def flag(value, good, warn, higher_is_better=True):
    """Return 'pass'/'warn'/'fail'/'na' for a value vs good/warn thresholds."""
    if value is None:
        return "na"
    if higher_is_better:
        return "pass" if value >= good else ("warn" if value >= warn else "fail")
    return "pass" if value <= good else ("warn" if value <= warn else "fail")


# Threshold flags (color-only; never fail the run). (good, warn, higher_is_better)
THRESHOLDS = {
    "alignment":       (90.0, 80.0, True),
    "mito":            (10.0, 20.0, False),
    "dup":             (20.0, 50.0, False),
    "FRiP":            (0.30, 0.20, True),
    "tss_enrichment":  (7.0, 5.0, True),
    "complexity_NRF":  (0.90, 0.80, True),
    "complexity_PBC1": (0.90, 0.80, True),
    "complexity_PBC2": (3.0, 1.0, True),
    "usable_fragments":(25_000_000, 15_000_000, True),   # ENCODE ATAC target >=25M
}


def _rd(path):
    try:
        with open(path) as fh:
            return fh.read()
    except OSError:
        return ""


def build_data(results_dir, samples):
    R = results_dir
    sections, summary = {}, {s: {} for s in samples}

    # Alignment (bowtie2, human index)
    rows = []
    raw_pairs = {}
    for s in samples:
        d = parse_bowtie2_log(_rd(f"{R}/aligned/{s}.bowtie2.log"))
        raw_pairs[s] = d.pop("total_pairs", None)   # raw read pairs; used by the usable-fragment section
        rows.append({"sample": s, **d})
        summary[s]["alignment"] = flag(d["overall_rate"], *THRESHOLDS["alignment"])
    sections["alignment"] = {"rows": rows, "flagkey": "overall_rate", "flagspec": "alignment",
                             "note": "vs the human index"}

    # Mitochondrial %
    rows = []
    for s in samples:
        pct = mito_pct_from_idxstats(_rd(f"{R}/filtered/{s}.idxstats.txt"))
        rows.append({"sample": s, "pct": pct})
        summary[s]["mito"] = flag(pct, *THRESHOLDS["mito"])
    sections["mito"] = {"rows": rows, "flagkey": "pct", "flagspec": "mito"}

    # Duplication (Picard)
    rows = []
    for s in samples:
        pct = parse_picard_dup(_rd(f"{R}/dedup/{s}.dedup.metrics.txt"))
        rows.append({"sample": s, "pct": pct})
        summary[s]["dup"] = flag(pct, *THRESHOLDS["dup"])
    sections["dup"] = {"rows": rows, "flagkey": "pct", "flagspec": "dup"}

    # Library complexity
    rows = []
    for s in samples:
        c = parse_complexity(_rd(f"{R}/library_complexity/{s}_complexity.txt"))
        rows.append({"sample": s, **c})
        for k in ("NRF", "PBC1", "PBC2"):
            summary[s][f"complexity_{k}"] = flag(c.get(k), *THRESHOLDS[f"complexity_{k}"])
    sections["complexity"] = {"rows": rows, "flagkey": "NRF", "flagspec": "complexity_NRF"}

    # Peaks + FRiP
    ps = read_tsv(f"{R}/qc/peak_summary.tsv") if os.path.exists(f"{R}/qc/peak_summary.tsv") else []
    for row in ps:
        try:
            summary[row["sample"]]["FRiP"] = flag(float(row["FRiP"]), *THRESHOLDS["FRiP"])
        except (KeyError, ValueError):
            pass
    sections["peaks"] = {"rows": ps, "flagkey": "FRiP", "flagspec": "FRiP"}

    # TSS enrichment
    tss = read_tsv(f"{R}/qc/tss_enrichment_scores.tsv") if os.path.exists(
        f"{R}/qc/tss_enrichment_scores.tsv") else []
    scorecol = "tss_enrichment" if (tss and "tss_enrichment" in tss[0]) else (
        list(tss[0].keys())[-1] if tss else None)
    for row in tss:
        if scorecol:
            try:
                summary[row["sample"]]["tss_enrichment"] = flag(
                    float(row[scorecol]), *THRESHOLDS["tss_enrichment"])
            except (KeyError, ValueError):
                pass
    sections["tss"] = {"rows": tss, "flagkey": scorecol, "flagspec": "tss_enrichment"}

    # Reads-in-annotation (tidy)
    ann = read_tsv(f"{R}/peak_annotation/reads_in_annotations.tsv") if os.path.exists(
        f"{R}/peak_annotation/reads_in_annotations.tsv") else []
    sections["annotation"] = {"rows": ann, "flagkey": None, "flagspec": None}

    # Blacklist filtering stats (whitespace-aligned text dump)
    sections["blacklist"] = {"rows": parse_blacklist_stats(_rd(f"{R}/qc/blacklist_filtering_stats.txt")),
                              "flagkey": None, "flagspec": None}

    # Usable fragments (final) + % of raw read pairs. The blacklist table's
    # Filtered_Reads is a samtools -c read count (both mates) -> fragments = reads/2;
    # raw read pairs come from the bowtie2 log's leading 'N reads' (PE = pairs).
    bl_by = {r["sample"]: r for r in sections["blacklist"]["rows"]}
    urows = []
    for s in samples:
        raw = raw_pairs.get(s)
        fr = bl_by.get(s, {}).get("filtered_reads")
        try:
            frag = int(fr) // 2
        except (TypeError, ValueError):
            frag = None
        pct = round(100.0 * frag / raw, 1) if (frag is not None and raw) else None
        urows.append({"sample": s, "raw_read_pairs": raw, "usable_fragments": frag, "usable_pct": pct})
        summary[s]["usable_fragments"] = flag(frag, *THRESHOLDS["usable_fragments"])
    sections["usable"] = {"rows": urows, "flagkey": "usable_fragments", "flagspec": "usable_fragments",
                          "note": "usable = unique, properly-paired, deduplicated, non-mito, non-blacklist "
                                  "fragments; usable_pct = usable fragments / raw read pairs. ENCODE ATAC target >=25M."}

    # Nucleosome signal from the fragment-length distribution (fragment_lengths.txt)
    fl = _rd(f"{R}/deeptools/fragment_lengths.txt")
    if fl.strip():
        by = {s["sample"]: s["points"] for s in parse_fragment_lengths(fl).get("series", [])}
        nrows = []
        for s in samples:
            st = nucleosome_stats(by.get(s, []))
            nrows.append({"sample": s, "nfr_pct": st["nfr_pct"], "mono_pct": st["mono_pct"],
                          "di_pct": st["di_pct"],
                          "mono_peak_bp": st["mono_peak_bp"] if st["mono_peak_bp"] is not None else "absent"})
            summary[s]["nfr_peak"] = st["nfr_peak"]
            summary[s]["mono_nucleosome"] = st["mono_nucleosome"]
        sections["nucleosome"] = {"rows": nrows, "flagkey": "mono_peak_bp", "flagspec": "mono_nucleosome",
            "note": "NFR < 147 bp; mononucleosome 147-294 bp (147..2x147); di 294-441 bp. ENCODE: a "
                    "nucleosome-free peak and a mononucleosome peak must both be present."}

    # Consensus counts (peak count)
    cbed = f"{R}/consensus/consensus_peaks.bed"
    n_consensus = sum(1 for ln in _rd(cbed).splitlines() if ln and not ln.startswith("#")) if os.path.exists(cbed) else 0
    sections["consensus"] = {"n_peaks": n_consensus}

    # Interactive charts (data-driven; no embedded images). Each guarded by file existence.
    colors = sample_colors(samples)
    D = f"{R}/deeptools"
    charts = {}

    def _attach(chart):
        for s in chart.get("series", []):
            s["color"] = colors.get(s["sample"], "#888888")
        return chart

    t = _rd(f"{D}/Profile_TSS.data.tab")
    if t.strip():
        charts["tss_profile"] = _attach(parse_tss_profile(t))

    t = _rd(f"{D}/ATACseq_fingerprint.tab")
    if t.strip():
        charts["fingerprint"] = _attach(parse_fingerprint(t))

    gc_series = []
    for s in samples:
        g = _rd(f"{D}/{s}.gc_content.txt")
        if g.strip():
            gc_series.append({"sample": s, "color": colors[s], "points": parse_gc_bias(g)})
    if gc_series:
        charts["gc_bias"] = {"x_label": "GC content (%)", "y_label": "Observed / expected",
                             "series": gc_series}

    t = _rd(f"{D}/deeptools_PCA.tab")
    if t.strip():
        pca = parse_pca(t)
        for p in pca["points"]:
            p["color"] = colors.get(p["sample"], "#888888")
        charts["pca"] = pca

    fl = _rd(f"{D}/fragment_lengths.txt")
    if fl.strip():
        charts["fragment_size"] = _attach(parse_fragment_lengths(fl))
    else:
        fs = _rd(f"{D}/fragmentsize.txt")
        if fs.strip():
            summ = parse_fragment_summary(fs)
            for row in summ["rows"]:
                row["color"] = colors.get(row["sample"], "#888888")
            charts["fragment_size"] = summ

    t = _rd(f"{D}/correlation_matrix.tab")
    if t.strip():
        cor = parse_cor_matrix(t)
        cor["colors"] = [colors.get(l, "#888888") for l in cor["labels"]]
        charts["correlation"] = cor

    t = _rd(f"{D}/tss_heatmap_downsampled.json")
    if t.strip():
        charts["tss_heatmap"] = parse_tss_heatmap(t)

    # Counting-unit label per section (kind drives the pill colour). "ratio" = the
    # metric is a same-unit ratio, so reads-vs-pairs cancels and the value is
    # unit-invariant; "frag" = molecule/fragment counting; "read" = per-read Tn5
    # insertions; "pairs" = read pairs (tool-native); "mixed" = per-column mix.
    UNITS = {
        "alignment":  ("mates / pairs (bowtie2)", "mixed"),
        "usable":     ("fragments", "frag"),
        "mito":       ("ratio · unit-free", "ratio"),
        "dup":        ("read pairs", "pairs"),
        "complexity": ("fragments", "frag"),
        "peaks":      ("reads · Tn5 cuts", "read"),
        "tss":        ("fragment coverage", "frag"),
        "nucleosome": ("fragments", "frag"),
        "annotation": ("reads · Tn5 cuts", "read"),
        "blacklist":  ("ratio · unit-free", "ratio"),
    }
    for k, (lab, kind) in UNITS.items():
        if k in sections:
            sections[k]["unit"] = lab
            sections[k]["unit_kind"] = kind

    return {"samples": samples, "thresholds": THRESHOLDS, "colors": colors,
            "sections": sections, "summary": summary, "charts": charts}


HEAD = """<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>ATAC-seq QC report</title>
<style>
:root{--bg:#fff;--fg:#1a1a1a;--muted:#666;--line:#e2e2e2;--surf:#f7f7f8;--bar:#7aa7d8;
--pass:#2e7d32;--warn:#b8860b;--fail:#c62828;--passbg:#e7f4e8;--warnbg:#fdf3dd;--failbg:#fbe6e6;}
@media (prefers-color-scheme:dark){:root{--bg:#15171a;--fg:#e8e8e8;--muted:#9aa0a6;--line:#2c2f33;
--surf:#1d2024;--bar:#4a6f9c;--passbg:#173a1c;--warnbg:#3a2f10;--failbg:#3a1717;}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);
font:14px/1.5 -apple-system,Segoe UI,Roboto,sans-serif;padding:24px;max-width:1100px;margin:auto}
h1{font-size:22px}h2{font-size:17px;border-bottom:1px solid var(--line);padding-bottom:4px;margin-top:32px}
table{border-collapse:collapse;width:100%;margin:8px 0;font-variant-numeric:tabular-nums}
th,td{padding:6px 10px;border-bottom:1px solid var(--line);text-align:right}
th:first-child,td:first-child{text-align:left}
th{cursor:pointer;user-select:none;background:var(--surf)}
.pass{color:var(--pass)}.warn{color:var(--warn)}.fail{color:var(--fail)}
.cell-pass{background:var(--passbg)}.cell-warn{background:var(--warnbg)}.cell-fail{background:var(--failbg)}
.bar{height:9px;background:var(--bar);border-radius:2px;display:inline-block;vertical-align:middle}
.muted{color:var(--muted)}img{max-width:100%;height:auto;border:1px solid var(--line);border-radius:4px;margin:6px 0}
details{margin:6px 0}summary{cursor:pointer;font-weight:600}
.unit{font-weight:500;font-size:11px;padding:1px 7px;border-radius:10px;margin-left:8px;
  vertical-align:middle;border:1px solid var(--line);white-space:nowrap}
.unit-frag{background:var(--passbg);color:var(--pass)}
.unit-read{background:#e6eef9;color:#1f5fae}
.unit-pairs{background:#efe7f7;color:#6a3fae}
.unit-ratio{background:var(--surf);color:var(--muted)}
.unit-mixed{background:var(--warnbg);color:var(--warn)}
@media (prefers-color-scheme:dark){.unit-read{background:#17273f;color:#8fb6ea}.unit-pairs{background:#241734;color:#b696e6}}
.unit-key{display:flex;flex-wrap:wrap;gap:8px;margin:8px 0 2px;font-size:12px}
.unit-key .unit{margin-left:0}
.charts{display:block}.chart{margin:14px 0}.chart svg{max-width:100%;height:auto}
.chart .ttl{font-weight:600;margin:10px 0 2px}
.legend{display:flex;flex-wrap:wrap;gap:10px;margin:6px 0}
.legend span{cursor:pointer;user-select:none;font-size:12px;display:inline-flex;align-items:center;gap:4px}
.legend span.off{opacity:.35;text-decoration:line-through}
.legend i{width:11px;height:11px;border-radius:2px;display:inline-block}
.axis{stroke:var(--line)}.axlab{fill:var(--muted);font-size:11px}
.tip{position:fixed;pointer-events:none;background:var(--surf);border:1px solid var(--line);
color:var(--fg);font-size:12px;padding:3px 6px;border-radius:4px;opacity:0;transition:opacity .1s;z-index:9}
canvas{max-width:100%;image-rendering:pixelated;border:1px solid var(--line);border-radius:4px}
</style></head><body>
<h1>ATAC-seq QC report</h1>
<p class="muted" id="meta"></p>
<div id="app"></div>
"""

JS = r"""
function el(t,c,h){var e=document.createElement(t);if(c)e.className=c;if(h!=null)e.innerHTML=h;return e;}
function short(s){return String(s).replace(/^GSF\d+-/,'');}
function fnum(x){return (x==null||x==='')?'':(typeof x==='number'?(Math.round(x*1000)/1000):x);}
function table(rows, cols, flagkey, flags){
  var t=el('table'); var thead=el('thead'); var tr=el('tr');
  cols.forEach(function(c){tr.appendChild(el('th',null,c));}); thead.appendChild(tr); t.appendChild(thead);
  var tb=el('tbody');
  var maxv={}; cols.forEach(function(c){rows.forEach(function(r){var v=parseFloat(r[c]);if(!isNaN(v))maxv[c]=Math.max(maxv[c]||0,v);});});
  rows.forEach(function(r){
    var row=el('tr');
    cols.forEach(function(c){
      var v=r[c]; var td=el('td');
      if(c===flagkey && flags && flags[r.sample]){td.className='cell-'+flags[r.sample];}
      var vf=parseFloat(v);
      if(c!=='sample' && !isNaN(vf) && maxv[c]){
        var w=Math.max(2,Math.round(100*vf/maxv[c]));
        td.innerHTML=fnum(v)+' <span class="bar" style="width:'+w+'px"></span>';
      } else { td.textContent = c==='sample'?short(v):fnum(v); }
      row.appendChild(td);
    });
    tb.appendChild(row);
  });
  t.appendChild(tb);
  thead.querySelectorAll('th').forEach(function(th,i){th.onclick=function(){sortBy(t,i);};});
  return t;
}
function sortBy(t,i){
  var tb=t.tBodies[0]; var rs=[].slice.call(tb.rows);
  var asc=t.getAttribute('data-sc')!=i+'a';
  rs.sort(function(a,b){var x=a.cells[i].textContent,y=b.cells[i].textContent;
    var nx=parseFloat(x),ny=parseFloat(y); if(!isNaN(nx)&&!isNaN(ny)){return asc?nx-ny:ny-nx;}
    return asc?x.localeCompare(y):y.localeCompare(x);});
  rs.forEach(function(r){tb.appendChild(r);}); t.setAttribute('data-sc',i+(asc?'a':'d'));
}
function unitPill(unit,kind){var p=el('span','unit unit-'+(kind||'ratio'),unit);
  p.title='Counting unit of this metric'+((kind==='ratio')?' — same-unit ratio, so reads vs pairs cancels (unit-free)':'');return p;}
function section(title, node, unit, kind){var d=el('details');d.open=true;
  var sm=el('summary',null,title); if(unit)sm.appendChild(unitPill(unit,kind));
  d.appendChild(sm); var body=el('div');body.appendChild(node);d.appendChild(body);return d;}
function colsOf(rows){return rows.length?Object.keys(rows[0]):[];}
var HIDDEN = {};                         // sample -> hidden? (shared across charts)
var TIP = null;
function tip(){ if(!TIP){TIP=el('div','tip');document.body.appendChild(TIP);} return TIP; }
function showTip(x,y,html){var t=tip();t.innerHTML=html;t.style.left=(x+12)+'px';t.style.top=(y+12)+'px';t.style.opacity=1;}
function hideTip(){if(TIP)TIP.style.opacity=0;}
var SVGNS='http:'+'//www.w3.org/2000/svg';   // split to avoid an external-looking URI substring in the static HTML; same value at runtime
function svgEl(t,a){var e=document.createElementNS(SVGNS,t);
  for(var k in a)e.setAttribute(k,a[k]);return e;}
function bounds(pts){var xs=pts.map(function(p){return p[0];}),ys=pts.map(function(p){return p[1];});
  return [Math.min.apply(null,xs),Math.max.apply(null,xs),Math.min.apply(null,ys),Math.max.apply(null,ys)];}

function frame(W,H,pad){
  var svg=svgEl('svg',{viewBox:'0 0 '+W+' '+H,width:W,height:H});
  return svg;
}
function axes(svg,W,H,pad,x0,x1,y0,y1,xLabel,yLabel){
  svg.appendChild(svgEl('line',{class:'axis',x1:pad.l,y1:H-pad.b,x2:W-pad.r,y2:H-pad.b}));
  svg.appendChild(svgEl('line',{class:'axis',x1:pad.l,y1:pad.t,x2:pad.l,y2:H-pad.b}));
  function tx(v){var t=svgEl('text',{class:'axlab','text-anchor':'middle',x:(pad.l+(W-pad.r))/2,y:H-4});t.textContent=xLabel;return t;}
  function ty(v){var t=svgEl('text',{class:'axlab','text-anchor':'middle',
    x:12,y:(pad.t+(H-pad.b))/2,transform:'rotate(-90 12 '+((pad.t+(H-pad.b))/2)+')'});t.textContent=yLabel;return t;}
  svg.appendChild(tx());svg.appendChild(ty());
  // numeric ticks (min/max)
  [[x0,pad.l],[x1,W-pad.r]].forEach(function(p){var t=svgEl('text',{class:'axlab','text-anchor':'middle',x:p[1],y:H-pad.b+14});t.textContent=(Math.round(p[0]*100)/100);svg.appendChild(t);});
  [[y0,H-pad.b],[y1,pad.t]].forEach(function(p){var t=svgEl('text',{class:'axlab','text-anchor':'end',x:pad.l-4,y:p[1]+3});t.textContent=(Math.round(p[0]*100)/100);svg.appendChild(t);});
}
function scaler(a0,a1,p0,p1){var d=(a1-a0)||1;return function(v){return p0+(v-a0)/d*(p1-p0);};}

function lineChart(cont,cfg){
  var W=680,H=300,pad={l:52,r:14,t:12,b:34};
  var vis=cfg.series.filter(function(s){return !HIDDEN[s.sample];});
  var all=[];vis.forEach(function(s){all=all.concat(s.points);});
  if(!all.length)all=[[0,0],[1,1]];
  var b=bounds(all);var sx=scaler(b[0],b[1],pad.l,W-pad.r),sy=scaler(b[2],b[3],H-pad.b,pad.t);
  var svg=frame(W,H,pad);axes(svg,W,H,pad,b[0],b[1],b[2],b[3],cfg.xLabel,cfg.yLabel);
  if(cfg.bands){cfg.bands.forEach(function(bd){
    if(bd.x<b[0]||bd.x>b[1])return;
    var xp=sx(bd.x);
    svg.appendChild(svgEl('line',{x1:xp.toFixed(1),y1:pad.t,x2:xp.toFixed(1),y2:H-pad.b,
      stroke:'#9993','stroke-width':1,'stroke-dasharray':'4 3'}));
    var lb=svgEl('text',{class:'axlab','text-anchor':'middle',x:xp.toFixed(1),y:pad.t+10,fill:'#999'});
    lb.textContent=bd.label;svg.appendChild(lb);
  });}
  cfg.series.forEach(function(s){
    if(HIDDEN[s.sample])return;
    var d=s.points.map(function(p,i){return (i?'L':'M')+sx(p[0]).toFixed(1)+' '+sy(p[1]).toFixed(1);}).join(' ');
    var path=svgEl('path',{d:d,fill:'none',stroke:s.color,'stroke-width':1.6});
    path.addEventListener('mousemove',function(e){
      showTip(e.clientX,e.clientY,'<b style="color:'+s.color+'">'+short(s.sample)+'</b>');});
    path.addEventListener('mouseleave',hideTip);
    svg.appendChild(path);
  });
  cont.appendChild(svg);
}
function scatterChart(cont,cfg){
  var W=520,H=380,pad={l:56,r:14,t:12,b:38};
  var pts=cfg.points.filter(function(p){return !HIDDEN[p.sample];});
  var xy=(pts.length?pts:cfg.points).map(function(p){return [p.x,p.y];});
  var b=bounds(xy);var mx=(b[1]-b[0])*0.08||1,my=(b[3]-b[2])*0.08||1;
  var sx=scaler(b[0]-mx,b[1]+mx,pad.l,W-pad.r),sy=scaler(b[2]-my,b[3]+my,H-pad.b,pad.t);
  var svg=frame(W,H,pad);axes(svg,W,H,pad,b[0]-mx,b[1]+mx,b[2]-my,b[3]+my,cfg.xLabel,cfg.yLabel);
  cfg.points.forEach(function(p){
    if(HIDDEN[p.sample])return;
    var c=svgEl('circle',{cx:sx(p.x).toFixed(1),cy:sy(p.y).toFixed(1),r:6,fill:p.color,stroke:'#0004'});
    c.addEventListener('mousemove',function(e){showTip(e.clientX,e.clientY,
      '<b style="color:'+p.color+'">'+short(p.sample)+'</b><br>'+cfg.xLabel+': '+fnum(p.x)+'<br>'+cfg.yLabel+': '+fnum(p.y));});
    c.addEventListener('mouseleave',hideTip);svg.appendChild(c);
    var tx=svgEl('text',{class:'axlab',x:(sx(p.x)+8).toFixed(1),y:(sy(p.y)+3).toFixed(1)});tx.textContent=short(p.sample);svg.appendChild(tx);
  });
  cont.appendChild(svg);
}
function heatmapGrid(cont,cfg){
  var n=cfg.labels.length,cell=Math.max(26,Math.min(60,360/n)),pad={l:96,t:96};
  var W=pad.l+n*cell+14,H=pad.t+n*cell+14;
  var svg=frame(W,H,{});
  function col(v){var t=(v+1)/2;var r=Math.round(255*(1-t)+33*t),g=Math.round(255*(1-t)+102*t),bl=Math.round(255*(1-t)+172*t);return 'rgb('+r+','+g+','+bl+')';}
  for(var i=0;i<n;i++){for(var j=0;j<n;j++){(function(i,j){
    var v=cfg.matrix[i][j];
    var rc=svgEl('rect',{x:pad.l+j*cell,y:pad.t+i*cell,width:cell-1,height:cell-1,fill:col(v)});
    rc.addEventListener('mousemove',function(e){showTip(e.clientX,e.clientY,
      short(cfg.labels[i])+' vs '+short(cfg.labels[j])+': <b>'+fnum(v)+'</b>');});
    rc.addEventListener('mouseleave',hideTip);svg.appendChild(rc);
  })(i,j);}}
  cfg.labels.forEach(function(l,k){
    var rt=svgEl('text',{class:'axlab','text-anchor':'end',x:pad.l-6,y:pad.t+k*cell+cell/2+3,fill:cfg.colors[k]});rt.textContent=short(l);svg.appendChild(rt);
    var ct=svgEl('text',{class:'axlab','text-anchor':'start',x:pad.l+k*cell+cell/2,y:pad.t-6,fill:cfg.colors[k],transform:'rotate(-45 '+(pad.l+k*cell+cell/2)+' '+(pad.t-6)+')'});ct.textContent=short(l);svg.appendChild(ct);
  });
  cont.appendChild(svg);
}
function canvasHeatmap(cont,cfg){
  var scale=2, gap=18;
  cfg.samples.forEach(function(s){
    var wrap=el('div');wrap.appendChild(el('div','axlab',short(s)));
    var cv=document.createElement('canvas');cv.width=cfg.ncols*scale;cv.height=cfg.nrows*scale;
    cv.style.width=Math.min(360,cfg.ncols*scale*2)+'px';
    var ctx=cv.getContext('2d');
    var bytes=atob(cfg.data[s]);var img=ctx.createImageData(cfg.ncols,cfg.nrows);
    for(var k=0;k<cfg.nrows*cfg.ncols;k++){var v=bytes.charCodeAt(k);var t=v/255;
      img.data[k*4]=Math.round(255*(1-t)+33*t);img.data[k*4+1]=Math.round(255*(1-t)+102*t);
      img.data[k*4+2]=Math.round(255*(1-t)+172*t);img.data[k*4+3]=255;}
    var tmp=document.createElement('canvas');tmp.width=cfg.ncols;tmp.height=cfg.nrows;
    tmp.getContext('2d').putImageData(img,0,0);ctx.imageSmoothingEnabled=false;
    ctx.drawImage(tmp,0,0,cfg.ncols*scale,cfg.nrows*scale);
    cv.title=short(s)+' — rows: '+cfg.nrows+' region-bins (sorted by signal), cols: TSS±'+cfg.upstream+'bp';
    wrap.appendChild(cv);wrap.style.display='inline-block';wrap.style.margin='0 '+gap+'px '+gap+'px 0';
    cont.appendChild(wrap);
  });
}
function legend(cont,samples,colors,onToggle){
  var lg=el('div','legend');
  samples.forEach(function(s){
    var sp=el('span',HIDDEN[s]?'off':null);
    sp.innerHTML='<i style="background:'+colors[s]+'"></i>'+short(s);
    sp.onclick=function(){HIDDEN[s]=!HIDDEN[s];sp.className=HIDDEN[s]?'off':'';onToggle();};
    lg.appendChild(sp);
  });
  cont.appendChild(lg);
}
function chartBlock(title,builder){
  var d=el('div','chart');d.appendChild(el('div','ttl',title));var host=el('div');d.appendChild(host);builder(host);return d;
}
function render(){
  document.getElementById('meta').textContent = 'Generated '+DATA.generated+' · '+DATA.samples.length+' samples';
  var app=document.getElementById('app');
  // Summary matrix
  var metrics=[]; Object.keys(DATA.summary).forEach(function(s){Object.keys(DATA.summary[s]).forEach(function(m){if(metrics.indexOf(m)<0)metrics.push(m);});});
  var st=el('table'); var hr=el('tr'); hr.appendChild(el('th',null,'sample'));
  metrics.forEach(function(m){hr.appendChild(el('th',null,m));}); var sh=el('thead');sh.appendChild(hr);st.appendChild(sh);
  var sb=el('tbody');
  DATA.samples.forEach(function(s){var tr=el('tr');tr.appendChild(el('td',null,short(s)));
    metrics.forEach(function(m){var f=(DATA.summary[s]||{})[m]||'na';var td=el('td','cell-'+f,f);tr.appendChild(td);});
    sb.appendChild(tr);}); st.appendChild(sb);
  app.appendChild(section('Summary (thresholds are color flags, not gates)', st));
  // Counting-unit key: what the tag next to each section title means
  var keyWrap=el('div');
  keyWrap.appendChild(el('p','muted','Each section is tagged with its counting unit — the unit of the denominator behind that metric:'));
  var key=el('div','unit-key');
  [['fragments','frag','molecule / fragment counting (dedup, complexity, insert size)'],
   ['reads · Tn5 cuts','read','per-read insertion counting (FRiP, promoter/enhancer)'],
   ['read pairs','pairs','tool-native pair counting (Picard duplication)'],
   ['ratio · unit-free','ratio','same-unit ratio: reads vs pairs cancels, value is unit-invariant'],
   ['mates / pairs (bowtie2)','mixed','bowtie2-native: overall rate is per-mate, concordant is per-pair']
  ].forEach(function(k){var w=el('span');w.appendChild(unitPill(k[0],k[1]));
    w.appendChild(document.createTextNode(' '+k[2]));key.appendChild(w);});
  keyWrap.appendChild(key);
  app.appendChild(section('Counting units (denominators)', keyWrap));
  // Numeric sections
  var order=[['alignment','Alignment rate'],['usable','Usable fragments (final)'],['mito','Mitochondrial %'],['dup','Duplication %'],
    ['complexity','Library complexity (NRF/PBC1/PBC2)'],
    ['peaks','Peaks + FRiP'],['tss','TSS enrichment'],['nucleosome','Nucleosome signal (fragment length)'],
    ['annotation','Reads in annotation'],['blacklist','Blacklist filtering']];
  order.forEach(function(o){var sec=DATA.sections[o[0]]; if(!sec||!sec.rows||!sec.rows.length)return;
    var flags={}; DATA.samples.forEach(function(s){var m=sec.flagspec; if(m&&DATA.summary[s]&&DATA.summary[s][m])flags[s]=DATA.summary[s][m];});
    var node=table(sec.rows, colsOf(sec.rows), sec.flagkey, flags);
    var wrap=el('div'); wrap.appendChild(node);
    if(sec.note)wrap.appendChild(el('p','muted',sec.note));
    app.appendChild(section(o[1], wrap, sec.unit, sec.unit_kind));
  });
  // Consensus
  if(DATA.sections.consensus){app.appendChild(section('Consensus peaks',
    el('p',null,'Consensus peak regions: <b>'+DATA.sections.consensus.n_peaks+'</b>')));}
  // Interactive charts (data-driven; consistent per-sample colors)
  var C=DATA.charts||{}, colors=DATA.colors||{};
  var cwrap=el('div','charts');
  var cnote=el('p','muted');cnote.appendChild(document.createTextNode('Coverage-based charts (TSS, fingerprint, correlation, PCA, GC) use '));
  cnote.appendChild(unitPill('fragment coverage','frag'));
  cnote.appendChild(document.createTextNode(' — bamCoverage --extendReads piles up whole fragments; the fragment-size / nucleosome chart counts '));
  cnote.appendChild(unitPill('fragments','frag'));cnote.appendChild(document.createTextNode('.'));
  cwrap.appendChild(cnote);
  legend(cwrap, DATA.samples, colors, redrawCharts);
  var host=el('div'); cwrap.appendChild(host);
  function redrawCharts(){
    host.innerHTML='';
    if(C.tss_profile) host.appendChild(chartBlock('TSS enrichment profile',function(h){lineChart(h,{series:C.tss_profile.series,xLabel:C.tss_profile.x_label,yLabel:C.tss_profile.y_label});}));
    if(C.fingerprint) host.appendChild(chartBlock('Fingerprint (cumulative signal)',function(h){lineChart(h,{series:C.fingerprint.series,xLabel:C.fingerprint.x_label,yLabel:C.fingerprint.y_label});}));
    if(C.fragment_size && C.fragment_size.mode==='dist') host.appendChild(chartBlock('Fragment-size distribution',function(h){lineChart(h,{series:C.fragment_size.series,xLabel:C.fragment_size.x_label,yLabel:C.fragment_size.y_label,bands:[{x:147,label:'NFR|mono'},{x:294,label:'mono|di'},{x:441,label:'di'}]});}));
    if(C.gc_bias) host.appendChild(chartBlock('GC bias',function(h){lineChart(h,{series:C.gc_bias.series,xLabel:C.gc_bias.x_label,yLabel:C.gc_bias.y_label});}));
    if(C.pca) host.appendChild(chartBlock('Sample PCA'+(C.pca.pc1_var!=null?' (PC1 '+C.pca.pc1_var+'% / PC2 '+C.pca.pc2_var+'%)':''),function(h){scatterChart(h,{points:C.pca.points,xLabel:'PC1',yLabel:'PC2'});}));
  }
  redrawCharts();
  // Non-toggle charts (correlation heatmap, canvas TSS heatmap, fragment summary fallback)
  if(C.correlation) cwrap.appendChild(chartBlock('Sample correlation (Spearman)',function(h){heatmapGrid(h,C.correlation);}));
  if(C.fragment_size && C.fragment_size.mode==='summary'){
    cwrap.appendChild(chartBlock('Fragment size (summary — re-run QC for full distribution)',function(h){
      var t=el('table'); var rows=C.fragment_size.rows;
      var head=el('tr');['sample','min','q1','median','q3','max'].forEach(function(c){head.appendChild(el('th',null,c));});
      var thead=el('thead');thead.appendChild(head);t.appendChild(thead);var tb=el('tbody');
      rows.forEach(function(r){var tr=el('tr');['sample','min','q1','median','q3','max'].forEach(function(c){tr.appendChild(el('td',null,c==='sample'?short(r[c]):fnum(r[c])));});tb.appendChild(tr);});
      t.appendChild(tb);h.appendChild(t);}));
  }
  if(C.tss_heatmap) cwrap.appendChild(chartBlock('TSS signal heatmap (per-region, downsampled)',function(h){canvasHeatmap(h,C.tss_heatmap);}));
  app.appendChild(section('Plots (interactive)', cwrap));
}
render();
"""


def render_html(data):
    return HEAD + "<script>\nconst DATA = " + json.dumps(data).replace("</", "<\\/") + ";\n" + JS + "\n</script>\n</body></html>\n"


def _discover_samples(results_dir):
    files = sorted(glob.glob(f"{results_dir}/library_complexity/*_complexity.txt"))
    return [os.path.basename(f)[: -len("_complexity.txt")] for f in files]


def main(argv=None):
    ap = argparse.ArgumentParser(description="Build the interactive ATAC-seq QC report.")
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--samples", default="", help="comma-separated; default = auto-discover")
    ap.add_argument("--generated", default="", help="timestamp string for the header")
    a = ap.parse_args(argv)
    samples = [s for s in a.samples.split(",") if s] or _discover_samples(a.results_dir)
    data = build_data(a.results_dir, samples)
    data["generated"] = a.generated
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write(render_html(data))
    print(f"[qc-report] wrote {a.out} ({len(samples)} samples, {len(data['charts'])} charts)")


if __name__ == "__main__":
    main()
