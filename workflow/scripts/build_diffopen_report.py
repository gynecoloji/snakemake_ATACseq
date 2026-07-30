#!/usr/bin/env python3
"""Self-contained HTML summary of the differential-openness run.

Renders every `results/diffopen/<mode>/` side by side, because the point of
running four normalizations is to see how much the answer depends on the choice.
No matplotlib / no external assets: charts are inline SVG.

Usage:
  python workflow/scripts/build_diffopen_report.py \\
      --diffopen-dir results/diffopen --out results/diffopen/diffopen_report.html
"""
from __future__ import annotations
import argparse, html, os
from datetime import datetime, timezone
import pandas as pd

CLASSES = [("all", "differential_openness"), ("promoter", "diffopen_promoter"),
           ("enhancer", "diffopen_enhancer")]
MODES = ("none", "ctcf", "rnastable")

CSS = """:root{--bg:#fff;--fg:#1a1a1a;--mut:#666;--line:#e3e3e3;--accent:#0b6e7c;
--warn:#b26a1b;--bad:#a12d2d;--good:#2d7a3e;--card:#fafafa}
@media(prefers-color-scheme:dark){:root{--bg:#16191c;--fg:#e8e8e8;--mut:#9aa0a6;--line:#2c3136;--card:#1d2125}}
*{box-sizing:border-box}body{margin:0;padding:2rem 1.25rem;background:var(--bg);color:var(--fg);
font:15px/1.6 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif}
.wrap{max-width:1060px;margin:0 auto}h1{font-size:1.6rem;margin:0 0 .25rem}
h2{font-size:1.15rem;margin:2.2rem 0 .6rem;padding-bottom:.3rem;border-bottom:2px solid var(--line)}
h3{font-size:1rem;margin:1.3rem 0 .4rem;color:var(--mut)}
.sub{color:var(--mut);font-size:.9rem;margin-bottom:1.5rem}
table{border-collapse:collapse;width:100%;margin:.6rem 0;font-size:.88rem}
th,td{padding:.42rem .6rem;text-align:right;border-bottom:1px solid var(--line)}
th:first-child,td:first-child{text-align:left}
th{font-weight:600;color:var(--mut);font-size:.78rem;text-transform:uppercase;letter-spacing:.04em}
tbody tr:hover{background:var(--card)}
code{background:var(--card);padding:.1rem .35rem;border-radius:3px;font-size:.85em}
.note{background:var(--card);border-left:3px solid var(--accent);padding:.8rem 1rem;margin:1rem 0;
border-radius:0 4px 4px 0;font-size:.9rem}.warn{border-left-color:var(--warn)}
.pill{display:inline-block;padding:.1rem .5rem;border-radius:10px;font-size:.78rem;font-weight:600}
.g{background:#2d7a3e22;color:var(--good)}.w{background:#b26a1b22;color:var(--warn)}
.b{background:#a12d2d22;color:var(--bad)}.scroll{overflow-x:auto}
footer{margin-top:3rem;padding-top:1rem;border-top:1px solid var(--line);color:var(--mut);font-size:.82rem}"""


def esc(x):
    return html.escape(str(x))


def read(p):
    return pd.read_csv(p, sep="\t") if os.path.exists(p) else None


def stats_of(d):
    if d is None or "pvalue" not in d.columns:
        return None
    nom, n01 = d[d.pvalue < 0.05], d[d.pvalue < 0.01]
    return dict(
        n=len(d),
        padj05=int((d.padj < 0.05).sum()) if "padj" in d.columns else 0,
        p05=len(nom), p01=len(n01),
        up05=100 * (nom.log2FoldChange > 0).mean() if len(nom) else float("nan"),
        up01=100 * (n01.log2FoldChange > 0).mean() if len(n01) else float("nan"))


def fmt(v, s="{:.1f}%"):
    return "&mdash;" if v is None or v != v else s.format(v)


def bars(labels, vals, ref=50.0, width=620, rowh=24):
    """Inline-SVG horizontal bars with a reference line (no plotting library)."""
    vmax = 100.0
    h = rowh * len(labels) + 24
    lw, tw = 200, width - 200 - 55
    o = [f'<svg viewBox="0 0 {width} {h}" width="100%" height="{h}" style="max-width:100%">']
    x = lw + tw * (ref / vmax)
    o.append(f'<line x1="{x:.0f}" y1="2" x2="{x:.0f}" y2="{h-18}" stroke="var(--mut)" '
             f'stroke-dasharray="4 3"/><text x="{x:.0f}" y="{h-5}" font-size="10" '
             f'fill="var(--mut)" text-anchor="middle">{ref:g}% (no bias)</text>')
    for i, (lab, v) in enumerate(zip(labels, vals)):
        y = i * rowh + 6
        o.append(f'<text x="0" y="{y+13}" font-size="11" fill="var(--fg)">{esc(lab)}</text>')
        if v != v:
            continue
        bw = max(tw * (v / vmax), 1)
        col = "var(--bad)" if abs(v - 100) < 1e-9 else "var(--accent)"
        o.append(f'<rect x="{lw}" y="{y+2}" width="{bw:.0f}" height="{rowh-9}" rx="2" '
                 f'fill="{col}" opacity=".85"/>'
                 f'<text x="{lw+bw+6:.0f}" y="{y+14}" font-size="11" fill="var(--mut)">{v:.1f}%</text>')
    return "".join(o) + "</svg>"


def top_table(d, n=8):
    if d is None or not len(d):
        return "<p class='sub'>none</p>"
    c = [x for x in ("Geneid", "Chr", "Start", "End", "log2FoldChange", "pvalue", "padj")
         if x in d.columns]
    s = d.nsmallest(n, "pvalue")[c]
    head = "".join(f"<th>{esc(x)}</th>" for x in c)
    rows = ""
    for _, r in s.iterrows():
        tds = "".join(
            f"<td>{r[x]:.3g}</td>" if isinstance(r[x], float) else f"<td>{esc(r[x])}</td>"
            for x in c)
        rows += f"<tr>{tds}</tr>"
    return f"<div class='scroll'><table><thead><tr>{head}</tr></thead><tbody>{rows}</tbody></table></div>"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--diffopen-dir", default="results/diffopen")
    ap.add_argument("--out", default="results/diffopen/diffopen_report.html")
    ap.add_argument("--contrast", default="NICD3 vs Control")
    a = ap.parse_args()
    root = a.diffopen_dir
    modes = [m for m in MODES if os.path.isdir(os.path.join(root, m))]

    st, sf = {}, {}
    for m in modes:
        st[m] = {c: stats_of(read(os.path.join(root, m, f"{s}.tsv"))) for c, s in CLASSES}
        t = read(os.path.join(root, m, "size_factors.tsv"))
        if t is not None:
            sf[m] = t.set_index("sample")["size_factor"]

    P = [f"<style>{CSS}</style><div class='wrap']".replace("]", "'>")]
    P.append("<h1>Differential openness &mdash; normalization comparison</h1>")
    P.append(f"<div class='sub'>{esc(a.contrast)} &middot; generated "
             f"{datetime.now(timezone.utc):%Y-%m-%d %H:%M UTC}</div>")

    # 1 - decision panel
    P.append("<h2>1. Which normalization should you trust?</h2><div class='scroll'><table>"
             "<thead><tr><th>mode</th><th>size-factor spread</th><th>peaks</th>"
             "<th>FDR&lt;0.05</th><th>p&lt;0.05</th><th>% up</th><th>verdict</th>"
             "</tr></thead><tbody>")
    for m in modes:
        s = st[m].get("all") or {}
        spr = (sf[m].max() / sf[m].min()) if m in sf else float("nan")
        up = s.get("up05", float("nan"))
        if up == up and abs(up - 100) < 1e-9:
            v = "<span class='pill b'>scaling artifact</span>"
        elif spr == spr and spr > 3:
            v = "<span class='pill w'>factors unstable</span>"
        elif spr != spr:
            v = "<span class='pill w'>n/a</span>"
        else:
            v = "<span class='pill g'>consistent</span>"
        P.append(f"<tr><td><code>{esc(m)}</code></td><td>{fmt(spr,'{:.2f}&times;')}</td>"
                 f"<td>{s.get('n',0):,}</td><td>{s.get('padj05','&mdash;')}</td>"
                 f"<td>{s.get('p05',0):,}</td><td>{fmt(up)}</td><td>{v}</td></tr>")
    P.append("</tbody></table></div>")
    P.append("<div class='note warn'><b>Read this first.</b> A large size-factor spread that "
             "also varies <i>within</i> a condition means the normalization is tracking handling "
             "variability, not biology. A <b>% up of exactly 100</b> is not a strong result &mdash; "
             "it is the fingerprint of a global scaling artifact. Prefer a mode where the "
             "internal methods (<code>none</code>, <code>ctcf</code>) agree.</div>")

    # 2 - size factors
    if sf:
        P.append("<h2>2. Per-sample size factors</h2><div class='scroll'><table><thead><tr>"
                 "<th>sample</th>" + "".join(f"<th>{esc(m)}</th>" for m in sf) +
                 "</tr></thead><tbody>")
        for smp in next(iter(sf.values())).index:
            P.append(f"<tr><td>{esc(smp)}</td>" +
                     "".join(f"<td>{sf[m].get(smp, float('nan')):.3f}</td>" for m in sf) + "</tr>")
        P.append("<tr><td><b>spread (max/min)</b></td>" +
                 "".join(f"<td><b>{sf[m].max()/sf[m].min():.2f}&times;</b></td>" for m in sf) +
                 "</tr></tbody></table></div>")

    # 3 - per class
    P.append("<h2>3. Results by peak class</h2>")
    P.append("<div class='note'>Promoter / enhancer assigned by overlap with the Ensembl "
             "Regulatory Build, <b>promoter precedence</b> (a peak hitting both counts once, as "
             "promoter). Each class is fit separately &rarr; own dispersion trend and within-class "
             "FDR; size factors are shared, computed once on all peaks.</div>")
    for m in modes:
        rows = [(c, st[m][c]) for c, _ in CLASSES if st[m].get(c)]
        if not rows:
            continue
        P.append(f"<h3>{esc(m)}</h3><div class='scroll'><table><thead><tr><th>class</th><th>n</th>"
                 "<th>padj&lt;0.05</th><th>p&lt;0.05</th><th>p&lt;0.01</th>"
                 "<th>% up (p&lt;0.05)</th><th>% up (p&lt;0.01)</th></tr></thead><tbody>")
        for c, s in rows:
            P.append(f"<tr><td>{esc(c)}</td><td>{s['n']:,}</td><td>{s['padj05']}</td>"
                     f"<td>{s['p05']:,}</td><td>{s['p01']:,}</td>"
                     f"<td>{fmt(s['up05'])}</td><td>{fmt(s['up01'])}</td></tr>")
        P.append("</tbody></table></div>")

    # 4 - direction balance
    P.append("<h2>4. Direction balance</h2>")
    labs, vals = [], []
    for m in modes:
        for c, _ in CLASSES:
            s = st[m].get(c)
            if s and s["p05"]:
                labs.append(f"{m} / {c}  (n={s['p05']:,})")
                vals.append(s["up05"])
    P.append(bars(labs, vals))
    P.append("<div class='note'><b>% up</b> = of the peaks reaching nominal p&lt;0.05, the fraction "
             "with log2FoldChange &gt; 0 (more open in treatment). It is a <b>direction balance, "
             "not a significance measure</b>. Near 50% = no coherent program. 70&ndash;95% = a real "
             "coordinated shift, which should <i>strengthen</i> as you tighten p&lt;0.05 &rarr; "
             "p&lt;0.01 (noise regresses to 50%, true signal sharpens). Exactly 100% = suspect a "
             "global scaling artifact.</div>")

    # 5 - top hits
    P.append("<h2>5. Strongest regions</h2>")
    for m in modes:
        for c, stem in CLASSES:
            d = read(os.path.join(root, m, f"{stem}.tsv"))
            if d is None or "pvalue" not in d.columns:
                continue
            if c == "all" and len(modes) > 1 and m != modes[-1]:
                pass
            P.append(f"<h3>{esc(m)} &middot; {esc(c)}</h3>{top_table(d)}")

    P.append("<footer>Generated by <code>workflow/scripts/build_diffopen_report.py</code>. "
             "With n=3 per condition the per-peak FDR is very conservative &mdash; treat "
             "nominal-p sets as a directional trend, not as individually defensible sites."
             "</footer></div>")

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write("<!doctype html><meta charset='utf-8'>"
                 "<meta name='viewport' content='width=device-width,initial-scale=1'>"
                 "<title>Differential openness report</title>" + "".join(P))
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
