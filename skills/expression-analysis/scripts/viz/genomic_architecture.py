"""Figure 1: genomic architecture of the regulatory window (footprint, Regulatory Build, transcript).

Every value comes from a core.locus.LocusConfig; nothing here is specific to one gene.
write_figure() writes figures/fig1_genomic_architecture.html: a header with the gene and
transcript, the figure as inline SVG (matplotlib, imported only when drawing), the footprint table
and the transcripts that start in the window.
"""

from __future__ import annotations

import html
import io
from pathlib import Path
from typing import Any

from core.locus import FootprintBlock, LocusConfig, footprint_blocks, load_transcripts_in_window
from viz.design_system import PALETTE, REGION_COLORS, mm_to_inches, region_color, region_label, set_publication_style

FIGURE_FILE = "figures/fig1_genomic_architecture.html"


def block_label(block: FootprintBlock) -> str:
    """Short label of a footprint block (two lines for exon parts)."""
    if block.kind == "5_prime_utr":
        return f"5'UTR\n({block.length}bp)"
    if block.kind == "cds":
        return f"CDS\n({block.length}bp)"
    return f"{block.kind.replace('_', ' ').capitalize()} ({block.length} bp)"


def region_name(block: FootprintBlock) -> str:
    return {
        "5_prime_utr": "Exon 1 (5' UTR)",
        "cds": "Exon 1 (Coding CDS)",
    }.get(block.kind, block.kind.replace("_", " ").capitalize())


def transcript_label(locus: LocusConfig) -> str:
    version = f".{locus.transcript_version}" if locus.transcript_version is not None else ""
    return f"{locus.transcript_id}{version}"


def primary_role(locus: LocusConfig, block: FootprintBlock) -> str:
    """Factual annotation source for the footprint table (no interpretive labels)."""
    if block.kind == "5_prime_utr":
        return f"{transcript_label(locus)} 5' UTR"
    if block.kind == "cds":
        five_prime = block.start if locus.strand == 1 else block.end
        first = abs(locus.distance_to_tss(five_prime) - locus.distance_to_tss(locus.cds_start_genomic))
        last = first + block.length - 1
        return f"{transcript_label(locus)} CDS codons {first // 3 + 1}-{last // 3 + 1}"
    if block.feature_id:
        return f"Ensembl Regulatory Build: {block.kind}"
    return "No Ensembl annotation"


def tick_positions(locus: LocusConfig, blocks: list[FootprintBlock]) -> list[int]:
    """Window start, every block start, CDS start (if inside the window) and window end."""
    ticks = {locus.window.start, locus.window.end, *(b.start for b in blocks)}
    if locus.window.contains(locus.cds_start_genomic):
        ticks.add(locus.cds_start_genomic)
    return sorted(ticks)


def labelled_ticks(locus: LocusConfig, ticks: list[int], min_spacing_fraction: float = 0.06) -> list[int]:
    """Subset of ticks that get a coordinate label, so that labels never overlap. Window ends are always labelled."""
    spacing = locus.window.length * min_spacing_fraction
    end = locus.window.end
    labelled = [ticks[0]]
    for pos in ticks[1:-1]:
        if pos - labelled[-1] >= spacing and end - pos >= spacing:
            labelled.append(pos)
    labelled.append(end)
    return labelled


def footprint_table_rows(locus: LocusConfig, blocks: list[FootprintBlock]) -> list[dict[str, str]]:
    chrom = locus.chromosome
    rows = [
        {
            "Region": region_name(b),
            "Ensembl ID": b.feature_id or "—",
            "Coordinates (GRCh38)": f"{chrom}:{b.start:,}-{b.end:,}",
            "Length": f"{b.length} bp",
            "SNVs (3/bp)": f"{b.snvs:,}",
            "Primary Role": primary_role(locus, b),
        }
        for b in blocks
    ]
    rows.append(
        {
            "Region": "TOTAL ISM UNIT",
            "Ensembl ID": "—",
            "Coordinates (GRCh38)": f"{chrom}:{locus.window.start:,}-{locus.window.end:,}",
            "Length": f"{locus.window.length} bp",
            "SNVs (3/bp)": f"{locus.window.length * 3:,}",
            "Primary Role": "Saturation mutagenesis window",
        }
    )
    return rows


def transcript_table_rows(locus: LocusConfig, transcripts: list[dict[str, str]]) -> list[dict[str, str]]:
    """Display rows of processed/transcripts_in_window.tsv."""
    rows = []
    for r in transcripts:
        version = f".{r['transcript_version']}" if r.get("transcript_version") else ""
        first = (f"{int(r['first_exon_start']):,}-{int(r['first_exon_end']):,}"
                 if r.get("first_exon_start") and r.get("first_exon_end") else "—")
        rows.append({
            "Transcript": f"{r['transcript_id']}{version}",
            "Biotype": r["transcript_type"],
            "MANE Select": "yes" if r["is_mane_select"] == "true" else "",
            "Selected": "yes" if r["is_selected"] == "true" else "",
            "TSS (GRCh38)": f"{locus.chromosome}:{int(r['tss']):,}",
            "Offset to selected TSS (bp)": r["tss_offset_bp"],
            "Exon 1": first,
        })
    return rows


def draw_regulatory_unit(locus: LocusConfig, blocks: list[FootprintBlock]):
    """Figure 1: genome-browser view of the regulatory unit (AlphaGenome / Nature figure style).

    Tracks share one GRCh38 axis: saturation-mutagenesis window, Ensembl Regulatory Build features,
    and the selected transcript (thin 5' UTR, thick CDS, intron backbone with strand chevrons,
    TSS arrow, start codon). Landmark coordinates label the axis; colours follow `design_system`.
    """
    set_publication_style()
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, Patch, Rectangle

    fig, ax = plt.subplots(figsize=mm_to_inches(180, 52))
    fig.subplots_adjust(left=0.15, right=0.98, top=0.96, bottom=0.27)
    win = locus.window
    margin = max(10, round(win.length * 0.02))
    lo, hi = win.start - margin, win.end + margin
    ax.set_xlim(lo, hi)
    axes_width_pt = fig.get_figwidth() * 72 * (0.98 - 0.15)
    bp_per_pt = (hi - lo) / axes_width_pt
    surface, muted, main = PALETTE["surface"], PALETTE["text_muted"], PALETTE["text_main"]

    def fits(text: str, width_bp: float, size: float) -> bool:
        return len(text) * size * 0.55 * bp_per_pt <= width_bp

    def ink_on(fill: str) -> str:
        r, g, b = (int(fill[i:i + 2], 16) / 255 for i in (1, 3, 5))
        lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
        return surface if lum < 0.5 else main

    def track_label(y: float, text: str) -> None:
        ax.text(-0.012, y, text, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=7, color=muted)

    ticks = tick_positions(locus, blocks)
    for pos in ticks:
        ax.axvline(pos, color=PALETTE["grid_line"], lw=0.5, zorder=0)

    # Regulatory Build lanes: overlapping features go to separate lanes
    features = sorted(locus.features, key=lambda f: (f.start, -f.end))
    lane_ends: list[int] = []
    lanes: list[int] = []
    for f in features:
        lane = next((i for i, end in enumerate(lane_ends) if end < f.start), len(lane_ends))
        if lane == len(lane_ends):
            lane_ends.append(f.end)
        lane_ends[lane] = f.end
        lanes.append(lane)
    n_lanes = max(1, len(lane_ends))
    y_reg = 1.55
    y_win = y_reg + 0.5 * (n_lanes - 1) + 0.85
    ax.set_ylim(-0.45, y_win + 0.55)

    # Track 1: saturation mutagenesis window (bracket)
    ax.plot([win.start, win.end], [y_win, y_win], color=PALETTE["ref"], lw=1.2, solid_capstyle="butt")
    for x in (win.start, win.end):
        ax.plot([x, x], [y_win - 0.12, y_win + 0.12], color=PALETTE["ref"], lw=1.2)
    ax.text((win.start + win.end) / 2, y_win + 0.18, f"{win.length:,} bp · {3 * win.length:,} SNVs",
            ha="center", va="bottom", fontsize=6.5, color=main)
    track_label(y_win, "Mutagenesis window")

    # Track 2: Ensembl Regulatory Build (labels inside the bar, when they fit)
    kinds: list[str] = []
    for f, lane in zip(features, lanes):
        s, e = max(f.start, lo), min(f.end, hi)
        y = y_reg + 0.5 * lane
        ax.add_patch(Rectangle((s - 0.5, y - 0.19), e - s + 1, 0.38, facecolor=region_color(f.biotype),
                               edgecolor=surface, lw=1.2, zorder=2))
        kinds.append(f.biotype)
        for text in (f"{region_label(f.biotype)} · {f.id}", region_label(f.biotype), f.id):
            if fits(text, e - s - 4 * bp_per_pt * 6.5, 6.5):
                ax.text((s + e) / 2, y, text, ha="center", va="center", fontsize=6.5,
                        color=ink_on(region_color(f.biotype)), zorder=3)
                break
    track_label(y_reg + 0.5 * (n_lanes - 1) / 2, "Regulatory Build")

    # Track 3: selected transcript (exon 1 parts inside the window, intron backbone towards 3')
    y_tx = 0.55
    ex1 = locus.first_exon
    three_prime = ex1.end if locus.strand == 1 else ex1.start
    edge = hi if locus.strand == 1 else lo
    ax.plot([three_prime, edge], [y_tx, y_tx], color=PALETTE["intron_backbone"], lw=0.8, zorder=1)
    marker = ">" if locus.strand == 1 else "<"
    step = max(1, min((hi - lo) / 45, abs(edge - three_prime) / 2))
    x = three_prime + step * locus.strand
    while (x < edge) if locus.strand == 1 else (x > edge):
        ax.plot(x, y_tx, marker=marker, ms=2.4, color=PALETTE["intron_backbone"], mew=0, zorder=1)
        x += step * locus.strand
    for b in blocks:
        if b.kind not in ("5_prime_utr", "cds"):
            continue
        h = 0.22 if b.kind == "5_prime_utr" else 0.44
        ax.add_patch(Rectangle((b.start - 0.5, y_tx - h / 2), b.length, h, facecolor=region_color(b.kind),
                               edgecolor=surface, lw=0.8, zorder=3))
        kinds.append(b.kind)
    tss = locus.tss
    y_arrow = y_tx + 0.62
    ax.plot([tss, tss], [y_tx + 0.11, y_arrow], color=PALETTE["tss"], lw=0.9, zorder=4)
    ax.add_patch(FancyArrowPatch((tss, y_arrow), (tss + 0.03 * (hi - lo) * locus.strand, y_arrow),
                                 arrowstyle="-|>", mutation_scale=6, color=PALETTE["tss"], lw=0.9, zorder=4))
    head = tss + 0.03 * (hi - lo) * locus.strand
    ax.text(head + 0.006 * (hi - lo) * locus.strand, y_arrow, f"TSS {tss:,}",
            ha="left" if locus.strand == 1 else "right", va="center", fontsize=6.5, color=main)
    if win.contains(locus.cds_start_genomic):
        sc_lo, sc_hi = locus.start_codon
        ax.add_patch(Rectangle((sc_lo - 0.5, y_tx - 0.22), sc_hi - sc_lo + 1, 0.44,
                               facecolor=region_color("start_codon"), edgecolor=surface, lw=0.6, zorder=4))
        kinds.append("start_codon")
        ax.text((sc_lo + sc_hi) / 2, y_tx - 0.32, "ATG", ha="center", va="top", fontsize=6.5, color=main)
    strand = "+" if locus.strand == 1 else "−"
    track_label(y_tx, f"{transcript_label(locus)} ({strand})")

    # Axis: landmark coordinates
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.set_yticks([])
    labelled = labelled_ticks(locus, ticks)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{p:,}" if p in labelled else "" for p in ticks])
    ax.set_xlabel(f"{locus.chromosome} position ({locus.assembly}, 1-based)")
    ax.spines["bottom"].set_position(("data", -0.3))

    seen = list(dict.fromkeys(k for k in REGION_COLORS if k in kinds))
    handles = [Patch(facecolor=region_color(k), edgecolor="none", label=region_label(k)) for k in seen]
    if len(handles) > 1:
        fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False,
                   bbox_to_anchor=(0.565, -0.02), handlelength=1.0, handleheight=0.8, columnspacing=1.4)
    return fig


def figure_svg(locus: LocusConfig, blocks: list[FootprintBlock]) -> str:
    """Figure 1 as an SVG document string (text kept as text, deterministic element IDs)."""
    import matplotlib.pyplot as plt

    fig = draw_regulatory_unit(locus, blocks)
    buffer = io.StringIO()
    try:
        fig.savefig(buffer, format="svg", metadata={"Date": None})
    finally:
        plt.close(fig)
    text = buffer.getvalue()
    return text[text.index("<svg"):]


def _table(rows: list[dict[str, Any]], caption: str) -> str:
    if not rows:
        return ""
    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in rows[0])
    body = "".join(
        "<tr>" + "".join(f"<td>{html.escape(str(v)).replace(chr(10), '<br>')}</td>" for v in row.values()) + "</tr>"
        for row in rows)
    return f"<table><caption>{html.escape(caption)}</caption><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


_CSS = f"""
body {{ font-family: Arial, Helvetica, sans-serif; color: {PALETTE['text_main']}; background: {PALETTE['surface']};
       margin: 0 auto; max-width: 1100px; padding: 16px; }}
h1 {{ font-size: 1.35rem; margin: 0 0 4px; }}
p.meta {{ color: {PALETTE['text_muted']}; margin: 2px 0; font-size: 0.9rem; }}
code {{ background: {PALETTE['bg_card']}; padding: 0 3px; border-radius: 3px; }}
figure {{ margin: 16px 0; overflow-x: auto; }}
figure svg {{ max-width: 100%; height: auto; }}
table {{ border-collapse: collapse; margin: 16px 0; font-size: 0.85rem; width: 100%; display: block; overflow-x: auto; }}
caption {{ text-align: left; font-weight: bold; padding: 4px 0; }}
th, td {{ border-bottom: 1px solid {PALETTE['border_frame']}; padding: 4px 8px; text-align: left; white-space: nowrap; }}
th {{ background: {PALETTE['bg_card']}; }}
"""


def figure_html(locus: LocusConfig, transcripts: list[dict[str, str]]) -> str:
    """Self-contained HTML page of Figure 1 with the footprint and transcript tables."""
    blocks = footprint_blocks(locus)
    e = html.escape
    strand = "forward strand" if locus.strand == 1 else "reverse strand"
    gene_version = f".{locus.gene_version}" if locus.gene_version is not None else ""
    decisions = "; ".join(f"{k} = {v}" for k, v in locus.decisions.items()) or "none"
    header = (
        f"<h1><i>{e(locus.gene_symbol)}</i> regulatory window and mutagenesis footprint</h1>"
        f"<p class='meta'>Gene <code>{e(locus.gene_symbol)}</code> (<code>{e(locus.gene_id + gene_version)}</code>) · "
        f"{e(locus.assembly)} {e(locus.chromosome)}, {strand}</p>"
        f"<p class='meta'>Transcript <code>{e(transcript_label(locus))}</code> ({e(locus.transcript_type)}"
        f"{', MANE Select' if locus.is_mane_select else ''}) · TSS <code>{e(locus.chromosome)}:{locus.tss:,}</code> · "
        f"CDS start <code>{e(locus.chromosome)}:{locus.cds_start_genomic:,}</code></p>"
        f"<p class='meta'>Window <code>{e(locus.window_str)}</code> ({locus.window.length:,} bp, "
        f"policy {e(locus.window_policy)}) · Decisions: {e(decisions)}</p>"
        "<p class='meta'>Gene structure: GENCODE v46 (alphagenome-variant-impact-score gtf, starts normalised to "
        "1-based). Regulatory features: Ensembl Regulatory Build (REST). Research use only.</p>"
    )
    tables = _table(footprint_table_rows(locus, blocks), "Regulatory unit footprint and coordinates")
    tables += _table(transcript_table_rows(locus, transcripts),
                     "Transcripts starting in the window (plus the selected transcript)")
    return (
        "<!doctype html>\n<html lang='en'><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width, initial-scale=1'>"
        f"<title>{e(locus.gene_symbol)} genomic architecture</title><style>{_CSS}</style></head><body>"
        f"{header}<figure>{figure_svg(locus, blocks)}</figure>{tables}</body></html>\n"
    )


def write_figure(study_dir: Path, locus: LocusConfig) -> Path:
    """Writes figures/fig1_genomic_architecture.html for the study; returns its path."""
    study_dir = Path(study_dir)
    path = study_dir / FIGURE_FILE
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(figure_html(locus, load_transcripts_in_window(study_dir)), encoding="utf-8")
    return path
