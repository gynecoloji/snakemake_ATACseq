#!/usr/bin/env python3
"""Downsample a deepTools computeMatrix (matrix.mat.gz) into a compact per-sample
byte heatmap JSON for the interactive QC report. Runs in the deepTools/numpy env."""
import argparse, base64, gzip, json
import numpy as np


def load_matrix(path):
    """matrix.mat.gz -> (header dict, 2D float array [regions, nsamp*bins]).
    First line is '@{json header}'; each data row has 6 BED cols then the values."""
    with gzip.open(path, "rt") as fh:
        header = json.loads(fh.readline()[1:])
        rows = []
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            rows.append([np.nan if x in ("nan", "") else float(x) for x in f[6:]])
    return header, np.array(rows, dtype=float)


def _bin_axis(a, n, axis):
    parts = [np.nanmean(np.take(a, ix, axis=axis), axis=axis)
             for ix in np.array_split(np.arange(a.shape[axis]), n)]
    return np.stack(parts, axis=axis)


def downsample(mat, sample_boundaries, nrows, ncols, vmax=None):
    """Sort regions by overall mean (desc, aligned across samples), average-bin each
    sample panel to nrows x ncols, quantize to uint8 on a shared vmax (~99th pct)."""
    order = np.argsort(-np.nan_to_num(np.nanmean(mat, axis=1), nan=-np.inf))
    mat = mat[order]
    if vmax is None:
        finite = mat[np.isfinite(mat)]
        vmax = float(np.nanpercentile(finite, 99)) if finite.size else 1.0
    if vmax <= 0:
        vmax = 1.0
    panels = {}
    for s in range(len(sample_boundaries) - 1):
        block = mat[:, sample_boundaries[s]:sample_boundaries[s + 1]]
        block = _bin_axis(_bin_axis(block, nrows, 0), ncols, 1)
        block = np.nan_to_num(block)
        panels[s] = (np.clip(block / vmax, 0, 1) * 255).astype(np.uint8)
    return panels, vmax


def main(argv=None):
    ap = argparse.ArgumentParser(description="Downsample computeMatrix to a byte heatmap JSON.")
    ap.add_argument("matrix")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--nrows", type=int, default=180)
    ap.add_argument("--ncols", type=int, default=80)
    a = ap.parse_args(argv)
    header, mat = load_matrix(a.matrix)
    panels, vmax = downsample(mat, header["sample_boundaries"], a.nrows, a.ncols)
    labels = header["sample_labels"]
    out = {"samples": labels, "nrows": a.nrows, "ncols": a.ncols,
           "upstream": int(header["upstream"][0]), "downstream": int(header["downstream"][0]),
           "vmax": vmax,
           "data": {labels[s]: base64.b64encode(panels[s].tobytes()).decode("ascii")
                    for s in range(len(labels))}}
    with open(a.out, "w") as fh:
        json.dump(out, fh)
    print(f"[downsample-tss] wrote {a.out} ({len(labels)} samples, {a.nrows}x{a.ncols})")


if __name__ == "__main__":
    main()
