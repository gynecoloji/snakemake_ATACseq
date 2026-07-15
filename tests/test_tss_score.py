import importlib.util, pathlib

_p = pathlib.Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "tss_score.py"
_spec = importlib.util.spec_from_file_location("tss_score", _p)
ts = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ts)


def test_enrichment_peaked_profile():
    # low background at the edges, high peak in the center
    vals = [1.0] * 10 + [10.0] + [1.0] * 10
    # center bin = 10, background mean = 1 -> enrichment ~10
    assert abs(ts.enrichment(vals) - 10.0) < 1e-6


def test_enrichment_flat_profile_is_one():
    assert abs(ts.enrichment([5.0] * 21) - 1.0) < 1e-6


def test_enrichment_empty_and_zero_background():
    assert ts.enrichment([]) == 0.0
    assert ts.enrichment([0.0, 0.0, 5.0, 0.0, 0.0]) == 0.0  # zero background -> 0


def test_parse_profile_skips_header_and_labels(tmp_path):
    f = tmp_path / "profile.tab"
    # header row of bin positions, then two sample rows with a group label column
    f.write_text(
        "bin\t-2000\t0\t2000\n"
        "SampleA\tgenes\t1.0\t8.0\t1.0\n"
        "SampleB\tgenes\t2.0\t2.0\t2.0\n"
    )
    prof = ts.parse_profile(str(f))
    assert set(prof) == {"SampleA", "SampleB"}
    assert prof["SampleA"] == [1.0, 8.0, 1.0]
    rows = ts.build(str(f))
    # SampleA is peaked (8 / 1), SampleB is flat (1.0)
    d = dict(rows)
    assert d["SampleA"] > d["SampleB"]
    assert abs(d["SampleB"] - 1.0) < 1e-6


def test_write_tsv(tmp_path):
    out = tmp_path / "tss.tsv"
    ts.write_tsv([("S1", 12.3), ("S2", 4.5)], out)
    lines = out.read_text().strip().splitlines()
    assert lines[0] == "sample\ttss_enrichment"
    assert lines[1] == "S1\t12.3"
