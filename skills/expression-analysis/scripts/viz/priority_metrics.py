"""Priority labels (colours and column documentation) for the workbench.

The labels themselves (`predicted_direction`, `target_min_abs_quantile`, `specificity_scope`) are
computed by core.contrast.add_priority_metrics and already stored in
processed/tissue_contrast.parquet; this module only names their values and documents the columns.
Nothing here is specific to one gene.
"""

from __future__ import annotations

from core.contrast import (  # noqa: F401  (re-exported for the workbench)
    DIRECTION_DOWN,
    DIRECTION_INCONSISTENT,
    DIRECTION_UP,
    SCOPE_MIXED,
    SCOPE_NO_CONTROL,
    SCOPE_NOT_CALLED,
    SCOPE_OPPOSITE,
    SCOPE_SHARED,
)

DIRECTION_COLORS = {
    DIRECTION_UP: "#16a34a",
    DIRECTION_DOWN: "#dc2626",
    DIRECTION_INCONSISTENT: "#f59e0b",
}


def column_docs(target_label: str = "target", control_label: str = "control") -> dict[str, str]:
    """Hover texts for workbench columns, with the cohort labels of the study."""
    t, c = target_label, control_label
    return {
        "variant_id": "GRCh38 SNV as chr:pos:ref>alt (1-based, forward strand).",
        "predicted_direction": f"Increase/Decrease when every {t} tissue's quantile in the direction contrast has that sign; otherwise Inconsistent.",
        "target_min_abs_quantile": f"Smallest |quantile| across {t} tissues of the direction contrast: the weakest {t} effect.",
        "specificity_scope": f"Votes of {c} tissues (one per tissue) against the {t} direction: Shared, Opposite or Mixed; No scope without a consistent direction or without {c} data.",
        "position": "1-based GRCh38 position.",
        "distance_to_tss": "Distance to the canonical TSS in transcript orientation (negative = upstream).",
        "genomic_domain": "Regulatory feature (upstream of the TSS) or Exon 1 part (5' UTR, start codon, CDS).",
        "ensembl_transcript_region": "Consequence class relative to the selected transcript.",
        "avi_phred": "AlphaGenome Variant Impact Phred: 10 = top 10%, 20 = top 1% of all SNVs.",
        "top_feature_key": "AVI feature with the largest SHAP attribution.",
        "target_mean": f"Mean raw effect over {t} tissues of a contrast (each tissue = mean of its plan tracks).",
        "control_mean": f"Mean raw effect over {c} tissues of a contrast.",
        "target_vs_control_diff": f"{t} mean minus {c} mean of a contrast's raw effects.",
        "top_tissue": "Tissue with the largest |raw effect| in the direction contrast, per cohort.",
        "atlas_url": "AlphaGenome Atlas page of the variant.",
        "ensembl_location_url": "Ensembl browser view around the variant.",
        "biosample": "Tissue name and ontology CURIE.",
        "cohort": f"{t} or {c} cohort of the study.",
        "assay": "Contrast (scorer and strand) a value comes from.",
        "raw_score": "AlphaGenome raw effect of ALT vs REF for this track.",
        "quantile_score": "Genome-wide quantile of the raw effect ([-1, 1] for signed scorers).",
    }


COLUMN_DOCS: dict[str, str] = column_docs()
