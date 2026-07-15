import importlib.util, pathlib, gzip, json, base64
import numpy as np
_p = pathlib.Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "downsample_tss_matrix.py"
_spec = importlib.util.spec_from_file_location("dts", _p)
dts = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(dts)


def _write_matrix(path):
    header = {"sample_labels": ["S1", "S2"], "sample_boundaries": [0, 4, 8],
              "group_boundaries": [0, 6], "upstream": [2000, 2000],
              "downstream": [2000, 2000], "bin size": [10, 10]}
    with gzip.open(path, "wt") as fh:
        fh.write("@" + json.dumps(header) + "\n")
        for r in range(6):                       # 6 regions, 6 BED cols + 8 value cols
            bed = ["chr1", str(r), str(r+1), f"g{r}", "0", "+"]
            vals = [str(float(r + c)) for c in range(8)]
            fh.write("\t".join(bed + vals) + "\n")


def test_load_and_downsample(tmp_path):
    m = tmp_path / "matrix.mat.gz"; _write_matrix(str(m))
    header, mat = dts.load_matrix(str(m))
    assert mat.shape == (6, 8)
    panels, vmax = dts.downsample(mat, header["sample_boundaries"], nrows=3, ncols=2)
    assert set(panels) == {0, 1}
    assert panels[0].shape == (3, 2) and panels[0].dtype == np.uint8
    assert vmax > 0


def test_cli_writes_json(tmp_path):
    m = tmp_path / "matrix.mat.gz"; _write_matrix(str(m))
    out = tmp_path / "hm.json"
    dts.main([str(m), "-o", str(out), "--nrows", "3", "--ncols", "2"])
    d = json.loads(out.read_text())
    assert d["samples"] == ["S1", "S2"]
    assert d["nrows"] == 3 and d["ncols"] == 2
    assert d["upstream"] == 2000 and d["downstream"] == 2000
    raw = base64.b64decode(d["data"]["S1"])
    assert len(raw) == 3 * 2            # nrows*ncols bytes
