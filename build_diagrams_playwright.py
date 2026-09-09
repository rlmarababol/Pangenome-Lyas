import json
import os
from playwright.sync_api import sync_playwright

with open(os.path.expanduser("~/Downloads/macrogen_1c/pangenome/pangenome_results/synteny/sequence_examples.json")) as f:
    data = json.load(f)

MATCH_COL = "#B491E5"
MISMATCH_COL = "#E4749B"
MATCH_TEXT = "#26215C"
MISMATCH_TEXT = "#4B1528"

CODON_STYLE_METRICS = {"frameshift", "inframe"}

METRIC_LABELS = {"snp": "SNP", "dN": "dN", "dS": "dS",
                  "frameshift": "Frameshift", "inframe": "In-frame indel"}

OUT_DIR = os.path.expanduser("~/Downloads/macrogen_1c/pangenome/pangenome_results/synteny/diagrams")
os.makedirs(OUT_DIR, exist_ok=True)

SKIP = {"inframe_shell_2"}  # 825bp indel, too large for a compact diagram


def html_wrapper(body, width_px):
    return f"""<html><head><style>
body {{ font-family: Helvetica, Arial, sans-serif; background: white; margin: 0; padding: 14px 16px; width: {width_px}px; }}
.title {{ font-weight: bold; font-size: 13px; margin-bottom: 2px; color: #1a1a1a; }}
.subtitle {{ font-size: 11px; color: #6b6b6b; margin-bottom: 12px; }}
.mono {{ font-family: 'Courier New', monospace; }}
.row {{ display: flex; align-items: center; line-height: 1.9; }}
.label {{ width: 130px; flex-shrink: 0; font-size: 11px; color: #6b6b6b; font-weight: bold; }}
.legend {{ display: flex; align-items: center; gap: 14px; margin-top: 12px; font-size: 11px; color: #6b6b6b; }}
.swatch {{ display: inline-block; width: 11px; height: 11px; border-radius: 2px; margin-right: 5px; vertical-align: middle; }}
</style></head><body>{body}</body></html>"""


def build_snp_html(key, d):
    ref_seq, qry_seq = d["ref_seq"], d["query_seq"]
    mismatches = set(d["mismatch_positions"])

    def render_seq(seq):
        spans = []
        for i, ch in enumerate(seq):
            if i in mismatches:
                spans.append(f'<span style="background:{MISMATCH_COL}66;color:{MISMATCH_TEXT};'
                             f'border-radius:2px;padding:0 1px;">{ch}</span>')
            else:
                spans.append(ch)
        return "".join(spans)

    body = f"""
    <div class="title">{d['family_id']} -- {METRIC_LABELS[key.split('_')[0]]} example</div>
    <div class="subtitle">{d['info']}</div>
    <div class="row"><span class="label">{d['ref_label']}</span><span class="mono">{render_seq(ref_seq)}</span></div>
    <div class="row"><span class="label">{d['query_label']}</span><span class="mono">{render_seq(qry_seq)}</span></div>
    <div class="legend"><span><span class="swatch" style="background:{MISMATCH_COL}66;"></span>SNP position</span></div>
    """
    width = max(500, len(ref_seq) * 9 + 180)
    return body, width


def build_codon_html(key, d):
    ref_nt, ref_aa = d["ref_nt"], d["ref_aa"]
    qry_nt, qry_aa = d["query_nt"], d["query_aa"]

    def codons(nt):
        return [nt[i:i+3] for i in range(0, len(nt) - len(nt) % 3, 3)]

    ref_codons = codons(ref_nt)
    qry_codons = codons(qry_nt)

    def render_row(aa_seq, codon_list, compare_aa):
        aa_cells, nt_cells = [], []
        for i, aa in enumerate(aa_seq):
            is_match = (compare_aa is None) or (i < len(compare_aa) and aa == compare_aa[i])
            bg = MATCH_COL if is_match else MISMATCH_COL
            fg = MATCH_TEXT if is_match else MISMATCH_TEXT
            aa_cells.append(f'<div style="min-width:22px;text-align:center;background:{bg}66;'
                             f'color:{fg};border-radius:3px;padding:2px 0;font-weight:bold;font-size:13px;">{aa}</div>')
            codon = codon_list[i] if i < len(codon_list) else ""
            nt_cells.append(f'<div style="min-width:22px;text-align:center;font-size:9px;color:#888;">{codon}</div>')
        return aa_cells, nt_cells

    ref_aa_cells, ref_nt_cells = render_row(ref_aa, ref_codons, None)
    qry_aa_cells, qry_nt_cells = render_row(qry_aa, qry_codons, ref_aa)

    def join_row(cells):
        return '<div class="mono" style="display:flex;gap:2px;">' + "".join(cells) + "</div>"

    body = f"""
    <div class="title">{d['family_id']} -- {METRIC_LABELS[key.split('_')[0]]} example</div>
    <div class="subtitle">{d['info']}</div>
    <div style="font-size:11px;color:#6b6b6b;font-weight:bold;margin-bottom:4px;">{d['ref_label']}</div>
    {join_row(ref_aa_cells)}
    {join_row(ref_nt_cells)}
    <div style="font-size:11px;color:#6b6b6b;font-weight:bold;margin:10px 0 4px;">{d['query_label']}</div>
    {join_row(qry_aa_cells)}
    {join_row(qry_nt_cells)}
    <div class="legend">
      <span><span class="swatch" style="background:{MATCH_COL}66;"></span>matches reference</span>
      <span><span class="swatch" style="background:{MISMATCH_COL}66;"></span>differs from reference</span>
    </div>
    """
    width = max(500, max(len(ref_aa), len(qry_aa)) * 24 + 60)
    return body, width


with sync_playwright() as p:
    browser = p.chromium.launch()
    # device_scale_factor renders at higher pixel density before screenshotting --
    # same idea as dpi=300 in every ggsave() call in this project. Default
    # screen rendering is ~96 DPI; scale_factor=3 gives ~288 DPI, close to the
    # 300 DPI standard used elsewhere. Bump to 4 for even crisper output.
    context = browser.new_context(device_scale_factor=3)
    page = context.new_page()

    n_saved = 0
    for key, d in data.items():
        if key in SKIP:
            print(f"Skipping {key} (too large for compact rendering)")
            continue
        metric = key.split("_")[0]
        if metric in CODON_STYLE_METRICS:
            body, width = build_codon_html(key, d)
        else:
            body, width = build_snp_html(key, d)

        png_path = f"{OUT_DIR}/{key}.png"
        page.set_viewport_size({"width": width, "height": 100})
        page.set_content(html_wrapper(body, width))
        height = page.evaluate("document.body.scrollHeight")
        page.set_viewport_size({"width": width, "height": height})
        page.screenshot(path=png_path)
        print(f"Rendered: {key}.png")
        n_saved += 1

    browser.close()

print(f"\nTotal rendered: {n_saved} / 17")
print(f"Saved in: {OUT_DIR}")
