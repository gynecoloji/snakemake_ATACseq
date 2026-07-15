import importlib.util, pathlib

_p = pathlib.Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "consensus_peaks.py"
_spec = importlib.util.spec_from_file_location("consensus_peaks", _p)
cp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cp)


def test_load_narrowpeak(tmp_path):
    f = tmp_path / "s_peaks.narrowPeak"
    f.write_text("chr1\t100\t300\tpeak1\t50\t.\t7.5\t9.1\t4.2\t80\n")
    peaks = cp.load_narrowpeak(str(f), "s")
    assert len(peaks) == 1
    p = peaks[0]
    assert (p.chrom, p.start, p.end, p.sample) == ("chr1", 100, 300, "s")
    assert abs(p.score - 4.2) < 1e-9
    assert p.summit == 180


def test_majority_keep_requires_two_replicates():
    P = cp.Peak
    rep1 = [P("chr1", 100, 300, 10, 200, "a"), P("chr1", 1000, 1100, 5, 1050, "a")]
    rep2 = [P("chr1", 150, 350, 10, 250, "b")]
    rep3 = []
    kept = cp.majority_keep([rep1, rep2, rep3], min_reps=2)
    summits = sorted(p.summit for p in kept)
    assert 200 in summits
    assert 250 in summits
    assert 1050 not in summits


def test_fixed_window_exact_width_and_clamp():
    P = cp.Peak
    assert cp.fixed_window(P("chr1", 900, 1100, 10, 1000, "a"), 500) == ("chr1", 750, 1250)
    c, s, e = cp.fixed_window(P("chr1", 0, 50, 10, 100, "a"), 500)
    assert (s, e) == (0, 500)


def test_assign_spm_is_per_sample():
    P = cp.Peak
    peaks = [P("chr1", 0, 1, 10, 0, "a"), P("chr1", 0, 1, 30, 0, "a"),
             P("chr1", 0, 1, 5, 0, "b")]
    cp.assign_spm(peaks)
    assert abs(peaks[0].spm - 250000.0) < 1e-3
    assert abs(peaks[2].spm - 1000000.0) < 1e-3


def test_iterative_overlap_removal_keeps_highest_spm():
    w = [
        {"chrom": "chr1", "start": 1000, "end": 1500, "spm": 100.0},
        {"chrom": "chr1", "start": 1200, "end": 1700, "spm": 50.0},
        {"chrom": "chr1", "start": 1600, "end": 2100, "spm": 40.0},
    ]
    kept = cp.iterative_overlap_removal(w, width=500)
    assert sorted(x["start"] for x in kept) == [1000, 1600]


def test_build_consensus_end_to_end(tmp_path):
    # 3-replicate majority group: one reproducible peak survives; a blacklisted
    # peak and a chrY peak and a single-replicate peak are all dropped.
    def wp(name, rows):
        p = tmp_path / name
        p.write_text("\n".join("\t".join(map(str, r)) for r in rows) + "\n")
        return p
    # narrowPeak cols: chr start end name score strand signal p q offset
    r1 = wp("A_peaks.narrowPeak", [
        ["chr1", 1000, 1400, "a1", 100, ".", 5, 5, 20.0, 200],   # summit 1200
        ["chr1", 5000, 5400, "a2", 100, ".", 5, 5, 10.0, 200],   # summit 5200 (blacklisted)
        ["chrY", 8000, 8400, "a3", 100, ".", 5, 5, 30.0, 200],   # dropped (chrY)
    ])
    r2 = wp("B_peaks.narrowPeak", [
        ["chr1", 1100, 1500, "b1", 100, ".", 5, 5, 18.0, 200],   # summit 1300 covers 1200
        ["chr1", 5050, 5450, "b2", 100, ".", 5, 5, 12.0, 150],   # summit 5200
    ])
    r3 = wp("C_peaks.narrowPeak", [
        ["chr1", 9000, 9400, "c1", 100, ".", 5, 5, 5.0, 200],    # only rep3 -> dropped
    ])
    blacklist = tmp_path / "bl.bed"
    blacklist.write_text("chr1\t5000\t5500\n")

    consensus = cp.build_consensus(
        {"g": ["A", "B", "C"]}, {"g": "majority"},
        {"A": str(r1), "B": str(r2), "C": str(r3)}, {},
        min_reps=2, width=500, keep_regex=r"^chr([1-9]|1[0-9]|2[0-2]|X)$",
        blacklist_path=str(blacklist),
    )
    assert len(consensus) == 1
    w = consensus[0]
    assert w["chrom"] == "chr1"
    assert w["end"] - w["start"] == 500
    assert w["start"] <= 1200 < w["end"]
    assert w["name"] == "consensus_peak_1"
