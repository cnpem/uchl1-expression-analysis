"""How names and numbers are written on the pages: genes, variants, coordinates and Portuguese numbers.

The rules are listed for the agent under *Writing conventions* in `references/glossary.md` (from
`assets/terms.yaml`); this module applies them where the code writes text, so every page writes them the same way.

* Human gene symbol: the HGNC symbol in italics (`*TP53*` in Markdown, `<i>TP53</i>` in HTML). The protein
  takes the same letters in roman type; the skill names the gene only.
* Variant: the Atlas identifier (`chr17:7676154:G>C`, the key every command and link uses) and, next to it, the
  genomic HGVS description on the RefSeq chromosome of GRCh38 (`NC_000017.11:g.7676154G>C`).
* Genomic coordinates: no thousands separator (`chr17:7676001-7676800`).
* Numbers in Portuguese text: decimal comma and a dot every three digits for counts (`1.234`, `12,5`).

Library: REFSEQ_GRCH38, gene, hgvs_g, count, dec.
"""

from __future__ import annotations

import re

# RefSeq accessions of the GRCh38 chromosomes (NCBI assembly report GCF_000001405.40, GRCh38.p14)
REFSEQ_GRCH38 = {
    "chr1": "NC_000001.11", "chr2": "NC_000002.12", "chr3": "NC_000003.12", "chr4": "NC_000004.12",
    "chr5": "NC_000005.10", "chr6": "NC_000006.12", "chr7": "NC_000007.14", "chr8": "NC_000008.11",
    "chr9": "NC_000009.12", "chr10": "NC_000010.11", "chr11": "NC_000011.10", "chr12": "NC_000012.12",
    "chr13": "NC_000013.11", "chr14": "NC_000014.9", "chr15": "NC_000015.10", "chr16": "NC_000016.10",
    "chr17": "NC_000017.11", "chr18": "NC_000018.10", "chr19": "NC_000019.10", "chr20": "NC_000020.11",
    "chr21": "NC_000021.9", "chr22": "NC_000022.11", "chrX": "NC_000023.11", "chrY": "NC_000024.10",
    "chrM": "NC_012920.1",
}
_SNV = re.compile(r"^(chr[0-9XYM]+):(\d+):([ACGT])>([ACGT])$")


def gene(symbol: str, fmt: str = "md") -> str:
    """The gene symbol in italics: Markdown (`md`), HTML (`html`) or as written (`plain`)."""
    if fmt == "html":
        return f"<i>{symbol}</i>"
    return f"*{symbol}*" if fmt == "md" else symbol


def hgvs_g(variant_id: str) -> str | None:
    """Genomic HGVS of an Atlas SNV identifier on GRCh38 (`chr17:7676154:G>C` -> `NC_000017.11:g.7676154G>C`);
    None for anything that is not a single-base substitution on a known chromosome."""
    m = _SNV.match(str(variant_id))
    if not m or m.group(1) not in REFSEQ_GRCH38:
        return None
    chrom, pos, ref, alt = m.groups()
    return f"{REFSEQ_GRCH38[chrom]}:g.{pos}{ref}>{alt}"


def count(n: int) -> str:
    """A count in Portuguese: a dot every three digits (2.382)."""
    return f"{int(n):,}".replace(",", ".")


def dec(x: float, sig: int = 3) -> str:
    """A value with `sig` significant digits and a decimal comma (21,4; 0,0121)."""
    return f"{x:.{sig}g}".replace(".", ",")
