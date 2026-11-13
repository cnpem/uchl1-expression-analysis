"""Shared visual tokens and publication styling for every figure (matplotlib and Altair).

Colour tokens follow the AlphaGenome figure conventions (REF slate vs ALT crimson, one hue per
assay modality) and Nature figure specs (89 / 180 mm widths). Region colours are one categorical
palette in a fixed order; the closest adjacent pair is labelled on the figures, so identity never
relies on colour alone. Nothing here is specific to one gene. matplotlib is imported only inside
set_publication_style.
"""

from __future__ import annotations

PALETTE: dict[str, str] = {
    # Allele contrast
    "ref": "#334155",
    "alt": "#dc2626",
    "ref_light": "#64748b",
    "alt_light": "#f87171",
    # Assay modalities
    "splice_sites": "#0891b2",
    "splice_site_usage": "#0284c7",
    "splice_junctions": "#4f46e5",
    "rna_seq": "#7c3aed",
    "cage": "#db2777",
    "dnase_atac": "#059669",
    "h3k27ac_enhancer": "#d97706",
    # Gene architecture
    "exon_cds": "#0369a1",
    "exon_utr": "#0ea5e9",
    "intron_backbone": "#94a3b8",
    "tss": "#dc2626",
    # Canvas
    "text_main": "#0f172a",
    "text_muted": "#475569",
    "grid_line": "#e2e8f0",
    "border_frame": "#cbd5e1",
    "bg_card": "#f8fafc",
    "surface": "#ffffff",
}

REGION_COLORS: dict[str, str] = {
    "enhancer": "#d97706",
    "promoter": "#4f46e5",
    "5_prime_utr": "#0ea5e9",
    "start_codon": "#db2777",
    "cds": "#0369a1",
    "open_chromatin_region": "#059669",
    "TF_binding_site": "#7c3aed",
    "CTCF_binding_site": "#0891b2",
}
"""Region class -> colour, in the fixed categorical order (genomic order of a typical regulatory unit)."""
OTHER_REGION_COLOR = "#94a3b8"

REGION_LABELS: dict[str, str] = {
    "enhancer": "Enhancer",
    "promoter": "Promoter",
    "5_prime_utr": "5' UTR",
    "start_codon": "Start codon",
    "cds": "CDS",
    "open_chromatin_region": "Open chromatin",
    "TF_binding_site": "TF binding site",
    "CTCF_binding_site": "CTCF binding site",
}


def region_color(kind: str) -> str:
    return REGION_COLORS.get(kind, OTHER_REGION_COLOR)


def region_label(kind: str) -> str:
    return REGION_LABELS.get(kind, kind.replace("_", " ").capitalize())


def domain_region_class(genomic_domain: str) -> str:
    """Region class of a `genomic_domain` label (core.locus.genomic_domain)."""
    text = genomic_domain.lower()
    if "start codon" in text:
        return "start_codon"
    if "5' utr" in text:
        return "5_prime_utr"
    if "cds" in text:
        return "cds"
    for kind in REGION_COLORS:
        if f"({kind.replace('_', ' ').lower()})" in text:
            return kind
    return "other"


def domain_color_scale(domains: list[str]) -> tuple[list[str], list[str]]:
    """Altair `domain` / `range` lists that colour `genomic_domain` values like Figure 1."""
    return list(domains), [region_color(domain_region_class(d)) for d in domains]


ASSAY_LABELS: dict[str, str] = {
    "RNA_SEQ": "RNA-seq",
    "RNA_SEQ_ACTIVE": "RNA-seq (active)",
    "CAGE": "CAGE",
    "PROCAP": "PRO-cap",
    "DNASE": "DNase-seq",
    "ATAC": "ATAC-seq",
    "CHIP_HISTONE": "histone ChIP-seq",
    "CHIP_TF": "TF ChIP-seq",
    "POLYADENYLATION": "polyadenylation",
    "SPLICE_SITES": "splice sites",
    "SPLICE_SITE_USAGE": "splice-site usage",
    "SPLICE_JUNCTIONS": "splice junctions",
}
"""AlphaGenome scorer -> assay name for axis titles and legends."""


def assay_label(scorer: str) -> str:
    return ASSAY_LABELS.get(scorer, scorer.replace("_", " ").lower())


# --------------------------------------------------------------------------- Altair (workbench figures)

PRINT_FONT = "Arial, Helvetica, sans-serif"


def apply_print_style(chart):
    """Top-level Altair config for paper / poster export: one sans font, larger text, recessive axes.

    Call it last, on the chart that is displayed (configure_* is only valid at the top level).
    """
    text, muted = PALETTE["text_main"], PALETTE["text_muted"]
    return (
        chart.configure(font=PRINT_FONT, background=PALETTE["surface"])
        .configure_axis(labelFontSize=13, titleFontSize=15, titleFontWeight="bold", labelColor=text, titleColor=text,
                        domainColor=muted, tickColor=muted, gridColor=PALETTE["grid_line"], labelLimit=320,
                        titlePadding=10)
        .configure_legend(labelFontSize=13, titleFontSize=14, titleColor=text, labelColor=text, symbolSize=140,
                          labelLimit=320, titleLimit=360, orient="bottom", direction="horizontal", columns=2,
                          columnPadding=18, rowPadding=6, titleOrient="top")
        .configure_title(fontSize=17, subtitleFontSize=13, color=text, subtitleColor=muted, anchor="start", offset=12)
        .configure_view(stroke=PALETTE["border_frame"])
    )


# --------------------------------------------------------------------------- matplotlib


def mm_to_inches(mm_width: float, mm_height: float) -> tuple[float, float]:
    return (mm_width / 25.4, mm_height / 25.4)


NATURE_SINGLE_COL = mm_to_inches(89, 70)
NATURE_INTERMEDIATE = mm_to_inches(120, 85)
NATURE_DOUBLE_COL = mm_to_inches(180, 105)

_RC_PARAMS = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7.0,
    "axes.labelsize": 7.0,
    "axes.titlesize": 8.0,
    "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5,
    "legend.title_fontsize": 7.0,
    "axes.linewidth": 0.6,
    "axes.edgecolor": PALETTE["text_main"],
    "axes.labelcolor": PALETTE["text_muted"],
    "xtick.color": PALETTE["text_muted"],
    "ytick.color": PALETTE["text_muted"],
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.top": False,
    "ytick.right": False,
    "xtick.minor.visible": False,
    "ytick.minor.visible": False,
    "xtick.major.size": 2.5,
    "xtick.major.width": 0.6,
    "ytick.major.size": 2.5,
    "ytick.major.width": 0.6,
    "grid.color": PALETTE["grid_line"],
    "grid.linestyle": "-",
    "grid.linewidth": 0.5,
    "figure.facecolor": PALETTE["surface"],
    "axes.facecolor": PALETTE["surface"],
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "svg.hashsalt": "fig1",
}


def set_publication_style() -> None:
    """SciencePlots `nature` + `no-latex` (when installed), then the shared rcParams on top.

    Uses the non-interactive Agg backend, so figures render without a display.
    """
    import matplotlib as mpl

    mpl.use("Agg")
    import matplotlib.pyplot as plt

    try:
        import scienceplots  # noqa: F401

        plt.style.use(["science", "nature", "no-latex"])
    except ImportError:
        pass
    mpl.rcParams.update(_RC_PARAMS)


__all__ = [
    "ASSAY_LABELS", "NATURE_DOUBLE_COL", "NATURE_INTERMEDIATE", "NATURE_SINGLE_COL", "OTHER_REGION_COLOR",
    "PALETTE", "REGION_COLORS", "REGION_LABELS", "apply_print_style", "assay_label", "domain_color_scale",
    "domain_region_class", "mm_to_inches", "region_color", "region_label", "set_publication_style",
]
