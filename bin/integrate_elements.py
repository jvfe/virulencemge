#!/usr/bin/env python3
"""Cross-reference virulence factors against prophages and genomic islands.

Takes the per-sample outputs of ABricate (VFDB), PhiSpy and IslandPath-DIMOB and
decides, for every virulence factor hit, whether it sits inside a predicted
pathogenicity island, inside a prophage, inside both, or in the core genome.

All coordinates are handled internally as 1-based inclusive intervals, matching
the convention used by ABricate, PhiSpy and GFF3. BED output is converted to
0-based half-open on the way out.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from urllib.parse import unquote

import pandas as pd

PROPHAGE = "prophage"
PAI = "pai"

CONTEXT_CORE = "Core genome"
CONTEXT_PAI = "PAI"
CONTEXT_PROPHAGE = "Prophage"
CONTEXT_BOTH = "PAI+Prophage"

VIRULENCE_COLUMNS = [
    "sample",
    "contig",
    "start",
    "end",
    "strand",
    "gene",
    "product",
    "database",
    "accession",
    "pct_coverage",
    "pct_identity",
    "gene_length",
    "mge_context",
    "in_pai",
    "in_prophage",
    "pai_ids",
    "prophage_ids",
    "pai_overlap_bp",
    "prophage_overlap_bp",
    "pai_overlap_fraction",
    "prophage_overlap_fraction",
]

ELEMENT_COLUMNS = [
    "sample",
    "element_id",
    "element_type",
    "contig",
    "start",
    "end",
    "length",
    "strand",
    "n_virulence_factors",
    "virulence_genes",
    "cross_element_ids",
    "cross_element_overlap_bp",
]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample", required=True, help="Sample identifier.")
    parser.add_argument("--abricate", type=Path, help="ABricate tabular report (VFDB).")
    parser.add_argument("--phispy", type=Path, help="PhiSpy prophage_coordinates.tsv.")
    parser.add_argument("--islandpath", type=Path, help="IslandPath-DIMOB GFF3 output.")
    parser.add_argument("--contig-map", type=Path, help="Contig identifier mapping table from prepare_genbank.py.")
    parser.add_argument(
        "--min-overlap",
        type=float,
        default=0.5,
        help="Minimum fraction of a virulence gene that must fall inside an element to call it contained (default: %(default)s).",
    )
    parser.add_argument("--prefix", help="Output file prefix (defaults to the sample identifier).")
    parser.add_argument("--outdir", type=Path, default=Path("."), help="Output directory.")
    return parser.parse_args(argv)


# --------------------------------------------------------------------------- #
# Contig identifier harmonisation
# --------------------------------------------------------------------------- #


def strip_version(identifier: str) -> str:
    return re.sub(r"\.\d+$", "", identifier)


class ContigResolver:
    """Translate tool-specific contig identifiers back to the FASTA identifier.

    PhiSpy reports the GenBank ACCESSION while IslandPath-DIMOB reports the
    LOCUS name, and either can differ from the FASTA header that ABricate
    reports. The mapping table written by ``prepare_genbank.py`` reconciles them;
    without it we fall back to comparing version-stripped identifiers.
    """

    def __init__(self, contig_map: Path | None = None) -> None:
        self.lookup: dict[str, str] = {}
        self.lengths: dict[str, int] = {}
        self.unresolved: set[str] = set()

        if contig_map and contig_map.exists():
            table = pd.read_csv(contig_map, sep="\t", dtype=str).fillna("")
            for _, row in table.iterrows():
                canonical = row["contig_id"]
                for column in ("contig_id", "locus_id", "accession_id"):
                    alias = row.get(column, "")
                    if alias:
                        self.lookup[alias] = canonical
                        self.lookup.setdefault(strip_version(alias), canonical)
                if row.get("length"):
                    self.lengths[canonical] = int(row["length"])

    def resolve(self, identifier: str) -> str:
        if identifier in self.lookup:
            return self.lookup[identifier]
        versionless = strip_version(identifier)
        if versionless in self.lookup:
            return self.lookup[versionless]
        if self.lookup:
            # A mapping table was supplied but this identifier is not in it.
            self.unresolved.add(identifier)
        return identifier


# --------------------------------------------------------------------------- #
# Input parsers
# --------------------------------------------------------------------------- #


def read_abricate(path: Path | None, sample: str, resolver: ContigResolver) -> pd.DataFrame:
    """Parse an ABricate report into normalised virulence factor hits."""
    columns = ["contig", "start", "end", "strand", "gene", "product", "database", "accession", "pct_coverage", "pct_identity"]
    if not path or not path.exists():
        return pd.DataFrame(columns=columns)

    table = pd.read_csv(path, sep="\t", dtype=str, comment=None).fillna("")
    table.columns = [column.lstrip("#").strip().upper() for column in table.columns]
    if "SEQUENCE" not in table.columns:
        return pd.DataFrame(columns=columns)

    hits = pd.DataFrame(
        {
            "contig": table["SEQUENCE"].map(resolver.resolve),
            "start": pd.to_numeric(table["START"], errors="coerce"),
            "end": pd.to_numeric(table["END"], errors="coerce"),
            "strand": table.get("STRAND", ""),
            "gene": table.get("GENE", ""),
            "product": table.get("PRODUCT", ""),
            "database": table.get("DATABASE", ""),
            "accession": table.get("ACCESSION", ""),
            "pct_coverage": pd.to_numeric(table.get("%COVERAGE"), errors="coerce"),
            "pct_identity": pd.to_numeric(table.get("%IDENTITY"), errors="coerce"),
        }
    )
    hits = hits.dropna(subset=["start", "end"])
    hits["start"] = hits["start"].astype(int)
    hits["end"] = hits["end"].astype(int)
    hits.insert(0, "sample", sample)
    return hits.reset_index(drop=True)


def read_phispy(path: Path | None, resolver: ContigResolver) -> list[dict]:
    """Parse the headerless PhiSpy prophage_coordinates.tsv.

    Columns are prophage number, contig, start and stop; att site columns may
    follow. Rows whose coordinates are not numeric are treated as a header.
    """
    elements: list[dict] = []
    if not path or not path.exists():
        return elements

    with open(path) as handle:
        for index, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            try:
                start, end = int(fields[2]), int(fields[3])
            except ValueError:
                continue
            elements.append(
                {
                    "element_id": fields[0] or f"pp{index}",
                    "element_type": PROPHAGE,
                    "contig": resolver.resolve(fields[1]),
                    "start": min(start, end),
                    "end": max(start, end),
                    "strand": ".",
                }
            )
    return elements


def read_islandpath(path: Path | None, resolver: ContigResolver) -> list[dict]:
    """Parse the GFF3 emitted by IslandPath-DIMOB."""
    elements: list[dict] = []
    if not path or not path.exists():
        return elements

    with open(path) as handle:
        for index, line in enumerate(handle, start=1):
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                continue

            attributes = {}
            if len(fields) > 8:
                for part in fields[8].split(";"):
                    key, _, value = part.partition("=")
                    if key:
                        attributes[key.strip()] = unquote(value.strip())

            elements.append(
                {
                    "element_id": attributes.get("ID") or f"GI_{index}",
                    "element_type": PAI,
                    "contig": resolver.resolve(fields[0]),
                    "start": min(start, end),
                    "end": max(start, end),
                    "strand": fields[6] if fields[6] in {"+", "-"} else ".",
                }
            )
    return elements


# --------------------------------------------------------------------------- #
# Interval arithmetic
# --------------------------------------------------------------------------- #


def overlap_bp(start_a: int, end_a: int, start_b: int, end_b: int) -> int:
    """Overlap of two 1-based inclusive intervals, in base pairs."""
    return max(0, min(end_a, end_b) - max(start_a, start_b) + 1)


def annotate_virulence_factors(
    hits: pd.DataFrame,
    elements: list[dict],
    min_overlap: float,
    sample: str,
) -> pd.DataFrame:
    """Assign a mobile genetic element context to every virulence factor hit."""
    by_contig: dict[str, list[dict]] = {}
    for element in elements:
        by_contig.setdefault(element["contig"], []).append(element)

    rows = []
    for hit in hits.to_dict("records"):
        gene_length = hit["end"] - hit["start"] + 1
        summary = {PAI: {"ids": [], "bp": 0}, PROPHAGE: {"ids": [], "bp": 0}}

        for element in by_contig.get(hit["contig"], []):
            shared = overlap_bp(hit["start"], hit["end"], element["start"], element["end"])
            if shared <= 0:
                continue
            bucket = summary[element["element_type"]]
            # Elements of one type can overlap each other, so track the largest
            # single overlap rather than summing double-counted base pairs.
            bucket["bp"] = max(bucket["bp"], shared)
            if shared / gene_length >= min_overlap:
                bucket["ids"].append(element["element_id"])

        in_pai = bool(summary[PAI]["ids"])
        in_prophage = bool(summary[PROPHAGE]["ids"])
        if in_pai and in_prophage:
            context = CONTEXT_BOTH
        elif in_pai:
            context = CONTEXT_PAI
        elif in_prophage:
            context = CONTEXT_PROPHAGE
        else:
            context = CONTEXT_CORE

        rows.append(
            {
                **hit,
                "sample": sample,
                "gene_length": gene_length,
                "mge_context": context,
                "in_pai": in_pai,
                "in_prophage": in_prophage,
                "pai_ids": ";".join(summary[PAI]["ids"]),
                "prophage_ids": ";".join(summary[PROPHAGE]["ids"]),
                "pai_overlap_bp": summary[PAI]["bp"],
                "prophage_overlap_bp": summary[PROPHAGE]["bp"],
                "pai_overlap_fraction": round(summary[PAI]["bp"] / gene_length, 4),
                "prophage_overlap_fraction": round(summary[PROPHAGE]["bp"] / gene_length, 4),
            }
        )

    return pd.DataFrame(rows, columns=VIRULENCE_COLUMNS)


def summarise_elements(
    elements: list[dict],
    annotated: pd.DataFrame,
    min_overlap: float,
    sample: str,
) -> pd.DataFrame:
    """Describe every predicted element, including the virulence genes it carries."""
    rows = []
    for element in elements:
        cross_type = PROPHAGE if element["element_type"] == PAI else PAI
        cross_ids, cross_bp = [], 0
        for other in elements:
            if other["element_type"] != cross_type or other["contig"] != element["contig"]:
                continue
            shared = overlap_bp(element["start"], element["end"], other["start"], other["end"])
            if shared > 0:
                cross_ids.append(other["element_id"])
                cross_bp += shared

        genes: list[str] = []
        if not annotated.empty:
            id_column = "pai_ids" if element["element_type"] == PAI else "prophage_ids"
            contained = annotated[
                (annotated["contig"] == element["contig"])
                & annotated[id_column].str.split(";").apply(lambda ids: element["element_id"] in ids)
            ]
            genes = sorted(set(contained["gene"]) - {""})

        rows.append(
            {
                "sample": sample,
                "element_id": element["element_id"],
                "element_type": element["element_type"],
                "contig": element["contig"],
                "start": element["start"],
                "end": element["end"],
                "length": element["end"] - element["start"] + 1,
                "strand": element["strand"],
                "n_virulence_factors": len(genes),
                "virulence_genes": ";".join(genes),
                "cross_element_ids": ";".join(cross_ids),
                "cross_element_overlap_bp": cross_bp,
            }
        )

    frame = pd.DataFrame(rows, columns=ELEMENT_COLUMNS)
    if frame.empty:
        return frame
    return frame.sort_values(["contig", "start", "element_type"]).reset_index(drop=True)


# --------------------------------------------------------------------------- #
# Output writers
# --------------------------------------------------------------------------- #

ELEMENT_LABELS = {PAI: "pathogenicity_island", PROPHAGE: "prophage"}


def write_element_tracks(elements: pd.DataFrame, resolver: ContigResolver, bed_path: Path, gff_path: Path) -> None:
    with open(bed_path, "w") as handle:
        handle.write('track name="mobile_genetic_elements" description="IslandPath-DIMOB PAIs and PhiSpy prophages"\n')
        for row in elements.to_dict("records"):
            name = f"{ELEMENT_LABELS[row['element_type']]}:{row['element_id']}"
            handle.write(f"{row['contig']}\t{row['start'] - 1}\t{row['end']}\t{name}\t0\t{row['strand']}\n")

    with open(gff_path, "w") as handle:
        handle.write("##gff-version 3\n")
        for contig in sorted(set(elements["contig"]) if not elements.empty else []):
            if contig in resolver.lengths:
                handle.write(f"##sequence-region {contig} 1 {resolver.lengths[contig]}\n")
        for row in elements.to_dict("records"):
            source = "islandpath" if row["element_type"] == PAI else "phispy"
            attributes = [
                f"ID={row['element_id']}",
                f"Name={ELEMENT_LABELS[row['element_type']]}_{row['element_id']}",
                f"element_type={ELEMENT_LABELS[row['element_type']]}",
                f"n_virulence_factors={row['n_virulence_factors']}",
            ]
            if row["virulence_genes"]:
                # GFF3 reserves ';' as the attribute separator; multiple values
                # within one attribute are comma-separated.
                attributes.append(f"virulence_genes={row['virulence_genes'].replace(';', ',')}")
            handle.write(
                "\t".join(
                    [
                        row["contig"],
                        source,
                        ELEMENT_LABELS[row["element_type"]],
                        str(row["start"]),
                        str(row["end"]),
                        ".",
                        row["strand"],
                        ".",
                        ";".join(attributes),
                    ]
                )
                + "\n"
            )


def write_virulence_track(annotated: pd.DataFrame, bed_path: Path) -> None:
    with open(bed_path, "w") as handle:
        handle.write('track name="virulence_factors" description="ABricate VFDB hits annotated with MGE context"\n')
        for row in annotated.to_dict("records"):
            strand = row["strand"] if row["strand"] in {"+", "-"} else "."
            name = f"{row['gene']}|{row['mge_context']}"
            handle.write(f"{row['contig']}\t{row['start'] - 1}\t{row['end']}\t{name}\t0\t{strand}\n")


def build_stats(sample: str, annotated: pd.DataFrame, elements: pd.DataFrame, resolver: ContigResolver) -> pd.DataFrame:
    contexts = annotated["mge_context"] if not annotated.empty else pd.Series(dtype=str)
    pais = elements[elements["element_type"] == PAI] if not elements.empty else elements
    prophages = elements[elements["element_type"] == PROPHAGE] if not elements.empty else elements

    n_vf = len(annotated)
    n_in_mge = int((contexts != CONTEXT_CORE).sum()) if n_vf else 0

    stats = {
        "sample": sample,
        "n_contigs": len(resolver.lengths) or annotated["contig"].nunique(),
        "n_virulence_factors": n_vf,
        "n_virulence_genes_unique": annotated["gene"].nunique() if n_vf else 0,
        "n_pai": len(pais),
        "pai_total_bp": int(pais["length"].sum()) if len(pais) else 0,
        "n_prophage": len(prophages),
        "prophage_total_bp": int(prophages["length"].sum()) if len(prophages) else 0,
        "n_vf_in_pai_only": int((contexts == CONTEXT_PAI).sum()) if n_vf else 0,
        "n_vf_in_prophage_only": int((contexts == CONTEXT_PROPHAGE).sum()) if n_vf else 0,
        "n_vf_in_pai_and_prophage": int((contexts == CONTEXT_BOTH).sum()) if n_vf else 0,
        "n_vf_in_core": int((contexts == CONTEXT_CORE).sum()) if n_vf else 0,
        "n_vf_in_mge": n_in_mge,
        "frac_vf_in_mge": round(n_in_mge / n_vf, 4) if n_vf else 0.0,
        "n_unresolved_contig_ids": len(resolver.unresolved),
    }
    return pd.DataFrame([stats])


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    prefix = args.prefix or args.sample
    args.outdir.mkdir(parents=True, exist_ok=True)

    resolver = ContigResolver(args.contig_map)

    hits = read_abricate(args.abricate, args.sample, resolver)
    elements = read_phispy(args.phispy, resolver) + read_islandpath(args.islandpath, resolver)

    annotated = annotate_virulence_factors(hits, elements, args.min_overlap, args.sample)
    element_table = summarise_elements(elements, annotated, args.min_overlap, args.sample)

    annotated.to_csv(args.outdir / f"{prefix}.virulence_annotation.tsv", sep="\t", index=False)
    element_table.to_csv(args.outdir / f"{prefix}.mge_regions.tsv", sep="\t", index=False)
    build_stats(args.sample, annotated, element_table, resolver).to_csv(
        args.outdir / f"{prefix}.stats.tsv", sep="\t", index=False
    )

    write_element_tracks(
        element_table,
        resolver,
        args.outdir / f"{prefix}.mge_regions.bed",
        args.outdir / f"{prefix}.mge_regions.gff3",
    )
    write_virulence_track(annotated, args.outdir / f"{prefix}.virulence_factors.bed")

    if resolver.unresolved:
        print(
            "WARNING: contig identifiers reported by PhiSpy or IslandPath-DIMOB that are absent from the "
            f"contig map: {', '.join(sorted(resolver.unresolved))}. Elements on these contigs cannot be "
            "cross-referenced against virulence factors.",
            file=sys.stderr,
        )

    print(
        f"{args.sample}: {len(annotated)} virulence factors, {len(element_table)} mobile elements, "
        f"{int((annotated['mge_context'] != CONTEXT_CORE).sum()) if len(annotated) else 0} virulence factors inside elements",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
