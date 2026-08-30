"""build -- inline each exported frame set into its page template.

    python -m analysis.model.first_principles.pages.build            # every page whose data exists
    python -m analysis.model.first_principles.pages.build ring cycle # named pages

Pages land in analysis/.pages/<name>.html (git-ignored; they embed megabytes
of frames). Export the frames first with the export_* modules, e.g.
`python -m analysis.model.first_principles.pages.export_ring`, and for the
overdrive `python -m analysis.model.first_principles.pages.export_overdrive
ring|hc15|block`. The pages are self-contained HTML: open them in a browser
or publish them as artifacts.
"""
from __future__ import annotations

import json
import sys

from analysis.model.first_principles.pages import common

PAGES = {
    "one_cell": ("one_cell", "one_cell"),
    "face_to_face": ("face_to_face", "face_to_face"),
    "vertex": ("vertex", "vertex"),
    "joints": ("joints", "joints"),
    "ring": ("ring", "ring"),
    "cycle": ("cycle", "cycle"),
    "overdrive_ring": ("overdrive", "overdrive_ring"),
    "overdrive_hc15": ("overdrive", "overdrive_hc15"),
    "overdrive_block": ("overdrive", "overdrive_block"),
}

# The overdrive template is written for the four-body ring; larger patches
# recolour by body kind and relabel the counts.
_PATCH = {
    "hc15": dict(title="Overdrive, Fifteen Bodies", cells=7, voids=8,
                 lede="The patch grown to one VE, its eight voids and its six axis VEs — fifteen bodies, thirty-two face welds — driven through the same 360°. Every body does what the ring's four did; at +60° all fifteen fold onto one octahedron.",
                 cross=(600, 612), ymax=650, scale=0.11, T=36000, scrub=120),
    "block": dict(title="Overdrive, Thirty-Five Bodies", cells=27, voids=8,
                  lede="The 3×3×3 block of VEs with its eight interior voids — thirty-five bodies, sixty-four face welds — driven through the same 360°. Every body does what the ring's four did; at +60° all thirty-five fold onto one octahedron.",
                  cross=(1272, 1620), ymax=1700, scale=0.075, T=40000, scrub=120),
}


def _patch_overdrive(t, patch):
    p = _PATCH[patch]
    nc, nv = p["cells"] * 8, p["voids"] * 8
    t = t.replace("<title>Overdrive</title>", f"<title>{p['title']}</title>")
    t = t.replace("<h1>Overdrive</h1>", f"<h1>{p['title']}</h1>")
    t = t.replace("The same ring, driven on past its dead end and kept going — a full 360° of the fold parameter.", p["lede"])
    t = t.replace("const COL = ['var(--o)','var(--u)','var(--d)','var(--l)'];",
                  "const COL = D.kind.map((k,i) => k==='cell' ? 'var(--o)' : 'var(--u)');")
    t = t.replace("isVoid:(k===1||k===3)", "isVoid:(D.kind[k]==='void')")
    t = t.replace('<span class="k">fronts facing out, cells</span><b class="mono" id="oc">8 · 8</b>\n      <span class="k">fronts facing out, voids</span><b class="mono" id="ov">8 · 8</b>',
                  f'<span class="k">fronts facing out, {p["cells"]} cells</span><b class="mono" id="oc">{nc} / {nc}</b>\n      <span class="k">fronts facing out, {p["voids"]} voids</span><b class="mono" id="ov">{nv} / {nv}</b>\n      <span class="k">spread of the centres</span><b class="mono" id="sp"></b>')
    t = t.replace("  document.getElementById('oc').textContent = fr.out[0] + ' · ' + fr.out[2];\n  document.getElementById('ov').textContent = fr.out[1] + ' · ' + fr.out[3];",
                  f"  document.getElementById('oc').textContent = fr.out_cells + ' / {nc}';\n  document.getElementById('ov').textContent = fr.out_voids + ' / {nv}';\n  document.getElementById('sp').textContent = fr.spread.toFixed(3);")
    t = t.replace('<span class="sw" style="background:var(--o)"></span><span>VE cell O · front</span>\n      <span class="sw" style="background:var(--d)"></span><span>VE cell D · front</span>\n      <span class="sw" style="background:var(--u)"></span><span>void U · front</span>\n      <span class="sw" style="background:var(--l)"></span><span>void L · front</span>',
                  f'<span class="sw" style="background:var(--o)"></span><span>the {p["cells"]} VEs · front</span>\n      <span class="sw" style="background:var(--u)"></span><span>the {p["voids"]} voids · front</span>')
    t = t.replace('<p class="s">Left scale: fronts out of 8, cells (indigo) and voids (teal). Right scale: crossing strut pairs among the four bodies (grey area) — 108 to 112 mid-passage, none at any multiple of 60°.</p>',
                  f'<p class="s">Left scale: fronts out as a fraction, cells (indigo) and voids (teal). Right scale: crossing strut pairs among all bodies (grey area) — {p["cross"][0]} to {p["cross"][1]} mid-passage, none at any multiple of 60°.</p>')
    t = t.replace("const yf = v => mt + (8 - v)/8*(H-mt-mb), yc = v => mt + (120 - v)/120*(H-mt-mb);",
                  f"const yf = v => mt + (1 - v)*(H-mt-mb), yc = v => mt + ({p['ymax']} - v)/{p['ymax']}*(H-mt-mb);")
    t = t.replace('for (const v of [0,4,8]) out += `<line x1="${ml}" x2="${W-mr}" y1="${yf(v)}" y2="${yf(v)}" stroke="var(--grid)"/><text x="${ml-8}" y="${yf(v)+4}" text-anchor="end" font-size="11" fill="var(--ink-2)" font-family="IBM Plex Mono,monospace">${v}</text>`;',
                  'for (const v of [0,0.5,1]) out += `<line x1="${ml}" x2="${W-mr}" y1="${yf(v)}" y2="${yf(v)}" stroke="var(--grid)"/><text x="${ml-8}" y="${yf(v)+4}" text-anchor="end" font-size="11" fill="var(--ink-2)" font-family="IBM Plex Mono,monospace">${v}</text>`;')
    t = t.replace("for (const v of [0,60,120]) out +=", f"for (const v of [0,{p['ymax']//2},{p['ymax']}]) out +=")
    t = t.replace("out += step(f=>f.out[0], 'var(--o)') + step(f=>f.out[1], 'var(--u)');",
                  f"out += step(f=>f.out_cells/{nc}, 'var(--o)') + step(f=>f.out_voids/{nv}, 'var(--u)');")
    t = t.replace("const S = Math.min(W,H) * 0.17 * zoom;", f"const S = Math.min(W,H) * {p['scale']} * zoom;")
    t = t.replace('max="180" step="1" value="0" aria-label="Fold of the VE cells, −60 to 300"', f'max="{p["scrub"]}" step="1" value="0" aria-label="Fold of the VE cells, −60 to 300"')
    t = t.replace("const T = 30000;", f"const T = {p['T']};")
    t = t.replace('stroke-width="1.3" stroke-linejoin="round"/>`;\n      else out += `<polygon', 'stroke-width="0.9" stroke-linejoin="round"/>`;\n      else out += `<polygon')
    return t


def build(name):
    template, data = PAGES[name]
    src = common.out(data)
    if not src.exists():
        print(f"  {name}: no frames at {src} -- run the exporter first")
        return False
    t = (common.TEMPLATES / f"{template}.html").read_text()
    if name.startswith("overdrive_") and name != "overdrive_ring":
        t = _patch_overdrive(t, name.split("_", 1)[1])
    html = t.replace("__DATA__", json.dumps(json.load(open(src)), separators=(",", ":")))
    assert "__DATA__" not in html
    common.PAGES.mkdir(parents=True, exist_ok=True)
    dst = common.PAGES / f"{name}.html"
    dst.write_text(html)
    print(f"  {name}: {len(html) // 1024} KB -> {dst}")
    return True


def main(argv):
    names = argv or list(PAGES)
    ok = [build(n) for n in names]
    return 0 if all(ok) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
