import importlib.util, pathlib
_p = pathlib.Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "build_qc_report.py"
_spec = importlib.util.spec_from_file_location("build_qc_report", _p)
bqr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bqr)

BOWTIE2 = """30000000 reads; of these:
  30000000 (100.00%) were paired; of these:
    3000000 (10.00%) aligned concordantly 0 times
    20000000 (66.67%) aligned concordantly exactly 1 time
    7000000 (23.33%) aligned concordantly >1 times
90.00% overall alignment rate
"""

IDXSTATS = "chr1\t248956422\t9000\t0\nchr2\t242193529\t8000\t0\nchrM\t16569\t3000\t0\n"

PICARD = """## htsjdk.samtools.metrics.StringHeader
# picard.sam.MarkDuplicates
LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE
Unknown Library\t0\t1000000\t0.0876\t5000000
"""

def test_parse_bowtie2_log():
    d = bqr.parse_bowtie2_log(BOWTIE2)
    assert d["overall_rate"] == 90.0
    assert d["concordant_uniq_pct"] == 66.67
    assert d["total_pairs"] == 30000000        # leading 'N reads' = raw read pairs (PE)

def test_mito_pct_from_idxstats():
    # 3000 chrM of 20000 total mapped = 15%
    assert round(bqr.mito_pct_from_idxstats(IDXSTATS), 2) == 15.0

def test_parse_picard_dup():
    assert round(bqr.parse_picard_dup(PICARD), 2) == 8.76

COMPLEXITY = """## Library Complexity Metrics for S1 ##
Total Reads\t100
NRF (Nd/Total)\t0.912000
PBC1 (N1/Nd)\t0.926000
PBC2 (N1/N2)\t15.26
"""

def test_parse_complexity():
    d = bqr.parse_complexity(COMPLEXITY)
    assert d["NRF"] == 0.912
    assert d["PBC1"] == 0.926
    assert d["PBC2"] == 15.26

def test_read_tsv(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("sample\tn_peaks\tFRiP\nS1\t1000\t0.31\nS2\t2000\t0.25\n")
    rows = bqr.read_tsv(str(p))
    assert rows[0] == {"sample": "S1", "n_peaks": "1000", "FRiP": "0.31"}
    assert rows[1]["n_peaks"] == "2000"

def test_flag():
    assert bqr.flag(0.95, 0.9, 0.8) == "pass"
    assert bqr.flag(0.85, 0.9, 0.8) == "warn"
    assert bqr.flag(0.7, 0.9, 0.8) == "fail"
    assert bqr.flag(5.0, 10.0, 20.0, higher_is_better=False) == "pass"   # mito %: lower better
    assert bqr.flag(25.0, 10.0, 20.0, higher_is_better=False) == "fail"
    assert bqr.flag(None, 0.9, 0.8) == "na"


BLACKLIST_STATS = """# Blacklist Filtering Statistics
# Date: Fri Jul 10 04:14:42 MDT 2026

                Sample  Original_Reads  Filtered_Reads  Blacklisted_Reads Percent_Excluded
 S1        54309380        54258408              50972            0.09%
 S2        76292828        76218820              74008            0.10%

Total reads before filtering: 130602208
Total reads after filtering: 130477228
Total blacklisted reads: 124980
Average percentage excluded: 0.10%
"""


def test_parse_blacklist_stats():
    rows = bqr.parse_blacklist_stats(BLACKLIST_STATS)
    assert len(rows) == 2
    assert rows[0] == {"sample": "S1", "original_reads": "54309380",
                        "filtered_reads": "54258408", "blacklisted_reads": "50972",
                        "pct_excluded": "0.09%"}
    assert rows[1]["sample"] == "S2"
    assert rows[1]["pct_excluded"] == "0.10%"
    assert all(r["sample"] != "Total" for r in rows)


def _mk_results(tmp_path):
    r = tmp_path / "results"
    (r / "aligned").mkdir(parents=True)
    (r / "filtered").mkdir(parents=True)
    (r / "dedup").mkdir(parents=True)
    (r / "library_complexity").mkdir(parents=True)
    (r / "qc").mkdir(parents=True)
    (r / "aligned" / "S1.bowtie2.log").write_text(BOWTIE2)
    (r / "filtered" / "S1.idxstats.txt").write_text(IDXSTATS)
    (r / "dedup" / "S1.dedup.metrics.txt").write_text(PICARD)
    (r / "library_complexity" / "S1_complexity.txt").write_text(COMPLEXITY)
    (r / "qc" / "peak_summary.tsv").write_text(
        "sample\tn_peaks\tmean_width\tmin_width\tmax_width\tFRiP\nS1\t100000\t480\t250\t900\t0.31\n")
    (r / "qc" / "tss_enrichment_scores.tsv").write_text("sample\ttss_enrichment\nS1\t8.2\n")
    (r / "qc" / "blacklist_filtering_stats.txt").write_text(BLACKLIST_STATS)
    (r / "deeptools").mkdir(parents=True)
    (r / "deeptools" / "Profile_TSS.data.tab").write_text(TSS_PROFILE)
    (r / "deeptools" / "S1.gc_content.txt").write_text(GC)
    (r / "deeptools" / "deeptools_PCA.tab").write_text(PCA)
    (r / "deeptools" / "ATACseq_fingerprint.tab").write_text(
        "#plotFingerprint --outRawCounts\n'a/S1.nobl.bam'\n" + "\n".join(str(i) for i in range(20)) + "\n")
    (r / "deeptools" / "correlation_matrix.tab").write_text(
        "\t'S1.nobl.bam'\n'S1.nobl.bam'\t1.0\n")
    (r / "deeptools" / "fragmentsize.txt").write_text(FRAG_TABLE)
    (r / "deeptools" / "fragment_lengths.txt").write_text(FRAG_LENGTHS)
    import json as _json, base64 as _b
    (r / "deeptools" / "tss_heatmap_downsampled.json").write_text(_json.dumps(
        {"samples": ["S1"], "nrows": 2, "ncols": 2, "upstream": 2000, "downstream": 2000,
         "vmax": 1.0, "data": {"S1": _b.b64encode(bytes([0, 64, 128, 255])).decode()}}))
    return str(r)

def test_build_data_shape(tmp_path):
    r = _mk_results(tmp_path)
    data = bqr.build_data(r, ["S1"])
    assert data["samples"] == ["S1"]
    assert data["sections"]["alignment"]["rows"][0]["overall_rate"] == 90.0
    assert round(data["sections"]["mito"]["rows"][0]["pct"], 1) == 15.0
    assert round(data["sections"]["dup"]["rows"][0]["pct"], 2) == 8.76
    assert data["sections"]["complexity"]["rows"][0]["NRF"] == 0.912
    assert data["sections"]["peaks"]["rows"][0]["FRiP"] == "0.31"
    # summary matrix carries a per-sample flag for each flagged metric
    assert data["summary"]["S1"]["complexity_NRF"] == "pass"
    assert data["summary"]["S1"]["mito"] == "warn"   # 15% is in 10-20 warn band
    assert data["sections"]["blacklist"]["rows"]
    # usable fragments: Filtered_Reads 54258408 / 2 = 27129204; raw pairs 30000000 -> 90.4%
    u = data["sections"]["usable"]["rows"][0]
    assert u["raw_read_pairs"] == 30000000
    assert u["usable_fragments"] == 27129204
    assert u["usable_pct"] == 90.4
    assert data["summary"]["S1"]["usable_fragments"] == "pass"    # >= 25M (ENCODE)
    # nucleosome signal: NFR peak + mononucleosome peak both present
    nuc = data["sections"]["nucleosome"]["rows"][0]
    assert nuc["mono_peak_bp"] == 200
    assert nuc["nfr_pct"] > nuc["mono_pct"] > nuc["di_pct"]
    assert data["summary"]["S1"]["nfr_peak"] == "pass"
    assert data["summary"]["S1"]["mono_nucleosome"] == "pass"
    # per-section counting-unit labels (denominator unit)
    assert data["sections"]["usable"]["unit_kind"] == "frag"
    assert data["sections"]["complexity"]["unit_kind"] == "frag"
    assert data["sections"]["nucleosome"]["unit_kind"] == "frag"
    assert data["sections"]["peaks"]["unit_kind"] == "read"        # FRiP = Tn5 insertions
    assert data["sections"]["annotation"]["unit_kind"] == "read"
    assert data["sections"]["dup"]["unit_kind"] == "pairs"
    assert data["sections"]["mito"]["unit_kind"] == "ratio"        # unit-invariant
    assert data["sections"]["blacklist"]["unit_kind"] == "ratio"
    assert data["sections"]["alignment"]["unit_kind"] == "mixed"
    assert data["sections"]["usable"]["unit"] == "fragments"
    assert data["colors"]["S1"] == "#4477aa"
    ch = data["charts"]
    assert ch["tss_profile"]["series"][0]["color"] == "#4477aa"
    assert ch["gc_bias"]["series"][0]["sample"] == "S1"
    assert ch["pca"]["points"][0]["sample"] == "S1"
    assert ch["fingerprint"]["series"][0]["points"]
    assert ch["fragment_size"]["mode"] in ("dist", "summary")
    assert ch["correlation"]["labels"] == ["S1"]
    assert ch["tss_heatmap"]["nrows"] == 2
    assert "plots" not in data                      # old embed path removed

def test_render_html_selfcontained(tmp_path):
    r = _mk_results(tmp_path)
    data = bqr.build_data(r, ["S1"])
    html = bqr.render_html(data)
    assert html.lstrip().startswith("<!DOCTYPE html>")
    assert "const DATA =" in html
    assert "http://" not in html and "https://" not in html   # no external refs
    assert "prefers-color-scheme" in html                     # theme-aware
    assert "S1" in html                                        # sample present in embedded DATA
    assert "data:image/png;base64," not in html     # nothing embedded as an image
    assert "<img" not in html                        # nothing linked as an image
    # per-metric counting-unit pills + legend key are rendered
    assert "unit-frag" in html and "unit-read" in html and "unit-ratio" in html
    assert "function unitPill" in html
    assert "Counting units (denominators)" in html

def test_main_writes_report(tmp_path):
    r = _mk_results(tmp_path)
    out = tmp_path / "report.html"
    bqr.main(["--results-dir", r, "--out", str(out), "--samples", "S1", "--generated", "2026-07-10"])
    html = out.read_text()
    assert html.startswith("<!DOCTYPE html>")
    assert "2026-07-10" in html and "const DATA =" in html

def test_sample_colors_stable_by_order():
    c = bqr.sample_colors(["Ctrl_1", "Ctrl_2", "NICD3_1"])
    assert c == {"Ctrl_1": "#4477aa", "Ctrl_2": "#66ccee", "NICD3_1": "#228833"}

def test_sample_colors_cycles_past_palette():
    many = [f"s{i}" for i in range(14)]
    c = bqr.sample_colors(many)
    assert c["s0"] == c["s12"]          # 12-color palette cycles
    assert len(set(c.values())) == 12

TSS_PROFILE = (
    "bin labels\t\t-2.0Kb\n"
    "bins\t\t1.0\t2.0\t3.0\t4.0\n"
    "S1\tgenes\t3.0\t3.2\t9.0\t3.1\n"
    "S2\tgenes\t2.0\t2.1\t6.0\t2.0\n"
)

def test_parse_tss_profile():
    d = bqr.parse_tss_profile(TSS_PROFILE)
    assert [s["sample"] for s in d["series"]] == ["S1", "S2"]
    s1 = d["series"][0]["points"]
    assert len(s1) == 4
    assert s1[0][1] == 3.0 and s1[2][1] == 9.0          # y values in order
    assert s1[0][0] < 0 < s1[-1][0]                     # x spans negative->positive around TSS

GC = ("0 100 1.0\n" "0 200 0.5\n" "0 300 1.5\n")       # 3 GC bins, ratio in col 3

def test_parse_gc_bias():
    pts = bqr.parse_gc_bias(GC)
    assert [p[1] for p in pts] == [1.0, 0.5, 1.5]
    assert pts[0][0] == 0.0 and pts[-1][0] == 100.0     # GC% spans 0..100

PCA = (
    "#plotPCA --outFileNameData\n"
    "Component\tS1.nobl.bam\tS2.nobl.bam\n"
    "1\t0.6\t-0.6\t75\n"
    "2\t0.4\t-0.4\t25\n"
)

def test_parse_pca():
    d = bqr.parse_pca(PCA)
    assert d["pc1_var"] == 75.0 and d["pc2_var"] == 25.0
    assert {p["sample"] for p in d["points"]} == {"S1", "S2"}
    p1 = next(p for p in d["points"] if p["sample"] == "S1")
    assert p1["x"] == 0.6 and p1["y"] == 0.4

PCA_REAL = (
    "#plotPCA --outFileNameData\n"
    "Component\tS1.nobl.bam\tS2.nobl.bam\tEigenvalue\n"
    "1\t0.6\t-0.6\t75\n"
    "2\t0.4\t-0.4\t25\n"
)

def test_parse_pca_with_eigenvalue_header():
    d = bqr.parse_pca(PCA_REAL)
    assert {p["sample"] for p in d["points"]} == {"S1", "S2"}   # no bogus 'Eigenvalue' point
    assert d["pc1_var"] == 75.0 and d["pc2_var"] == 25.0

def test_parse_fingerprint_cumulative_downsampled():
    header = "'a/S1.nobl.bam'\t'a/S2.nobl.bam'"
    body = "\n".join(f"{i}\t{2*i}" for i in range(1000))     # S2 = 2*S1, both increasing
    txt = "#plotFingerprint --outRawCounts\n" + header + "\n" + body + "\n"
    d = bqr.parse_fingerprint(txt, n_points=50)
    assert [s["sample"] for s in d["series"]] == ["S1", "S2"]
    s1 = d["series"][0]["points"]
    assert len(s1) <= 50
    assert s1[0][0] == 0.0 and abs(s1[-1][0] - 1.0) < 1e-6      # x = rank fraction 0..1
    assert abs(s1[-1][1] - 1.0) < 1e-6                          # y = cumulative fraction ->1
    assert all(s1[i][1] <= s1[i+1][1] + 1e-9 for i in range(len(s1)-1))  # monotone

FRAG_RAW = ("#bamPEFragmentSize\n"
            "Size\tOccurrences\tSample\n"
            "50\t5\tS1.nobl.bam\n"
            "60\t9\tS1.nobl.bam\n"
            "50\t3\tS2.nobl.bam\n")

def test_parse_fragment_lengths():
    d = bqr.parse_fragment_lengths(FRAG_RAW)
    assert d["mode"] == "dist"
    assert [s["sample"] for s in d["series"]] == ["S1", "S2"]
    assert d["series"][0]["points"] == [[50.0, 5.0], [60.0, 9.0]]


def _frag_lengths_text(*samples):
    # Synthetic distribution with an NFR peak at ~60 bp and a mononucleosome peak
    # at ~200 bp, plus a small flat baseline — a "good" ATAC fragment profile.
    lines = ["#bamPEFragmentSize", "Size\tOccurrences\tSample"]
    for s in samples:
        for L in range(20, 400):
            nfr = max(0, 100 - abs(L - 60) * 2)      # sub-nucleosomal peak at 60 bp
            mono = max(0, 50 - abs(L - 200))         # mononucleosome peak at 200 bp
            lines.append(f"{L}\t{nfr + mono + 2}\t{s}.nobl.bam")
    return "\n".join(lines) + "\n"

FRAG_LENGTHS = _frag_lengths_text("S1")

def test_nucleosome_stats_detects_both_peaks():
    d = bqr.parse_fragment_lengths(FRAG_LENGTHS)
    st = bqr.nucleosome_stats(d["series"][0]["points"])
    assert st["nfr_peak"] == "pass"                  # global mode is sub-nucleosomal (< 147 bp)
    assert st["mono_nucleosome"] == "pass"           # a mononucleosome bump is present
    assert st["mono_peak_bp"] == 200                 # detected at the 200-bp peak
    assert st["nfr_pct"] > st["mono_pct"] > st["di_pct"]
    assert abs(st["nfr_pct"] + st["mono_pct"] + st["di_pct"] - 100.0) < 1.5

def test_nucleosome_stats_empty():
    st = bqr.nucleosome_stats([])
    assert st == {"nfr_pct": None, "mono_pct": None, "di_pct": None, "mono_peak_bp": None,
                  "nfr_peak": "na", "mono_nucleosome": "na"}

def test_nucleosome_stats_no_mono_peak():
    # Monotonically decaying distribution: NFR only, no mononucleosome bump.
    pts = [[float(L), float(max(1, 300 - L))] for L in range(20, 300)]
    st = bqr.nucleosome_stats(pts)
    assert st["nfr_peak"] == "pass"
    assert st["mono_nucleosome"] == "fail"
    assert st["mono_peak_bp"] is None

FRAG_TABLE = (
    "\tFrag. Sampled\tFrag. Len. Min.\tFrag. Len. 1st. Qu.\tFrag. Len. Mean\t"
    "Frag. Len. Median\tFrag. Len. 3rd Qu.\tFrag. Len. Max\n"
    "a/S1.nobl.bam\t28467\t31.0\t55.0\t110.3\t81.0\t133.0\t964.0\n")

def test_parse_fragment_summary():
    d = bqr.parse_fragment_summary(FRAG_TABLE)
    assert d["mode"] == "summary"
    r = d["rows"][0]
    assert r["sample"] == "S1"
    assert r["median"] == 81.0 and r["q1"] == 55.0 and r["q3"] == 133.0
    assert r["min"] == 31.0 and r["max"] == 964.0

CORMAT = ("#plotCorrelation --outFileCorMatrix\n"
          "\t'S1.nobl.bam'\t'S2.nobl.bam'\n"
          "'S1.nobl.bam'\t1.0\t0.8\n"
          "'S2.nobl.bam'\t0.8\t1.0\n")

def test_parse_cor_matrix():
    d = bqr.parse_cor_matrix(CORMAT)
    assert d["labels"] == ["S1", "S2"]
    assert d["matrix"] == [[1.0, 0.8], [0.8, 1.0]]

def test_html_has_chart_primitives(tmp_path):
    r = _mk_results(tmp_path)
    html = bqr.render_html(bqr.build_data(r, ["S1"]))
    for fn in ["function lineChart", "function scatterChart", "function heatmapGrid",
               "function canvasHeatmap", "function legend", "DATA.charts"]:
        assert fn in html
    assert "document.createElement('canvas')" in html   # heatmap uses canvas
    assert "data:image/png;base64," not in html
