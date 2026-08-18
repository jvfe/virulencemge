#!/usr/bin/env python3
"""Build a canonical GenBank file for prophage and genomic island prediction.

PhiSpy and IslandPath-DIMOB both consume GenBank, but they report contig
identifiers differently: PhiSpy uses the record ``id`` (the ACCESSION/VERSION
line) while IslandPath-DIMOB uses the record ``name`` (the LOCUS line). ABricate
in turn reports the FASTA header. This script emits a GenBank file in which all
three are derived from the same FASTA contig identifier, plus a mapping table
that lets downstream code recover the FASTA identifier when the LOCUS line had
to be shortened to satisfy the GenBank format.

It accepts either an existing GenBank file (which is harmonised) or a GFF3
annotation (from which a GenBank file is built), always adding the nucleotide
sequence and CDS translations that IslandPath-DIMOB requires.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
import sys
from collections import OrderedDict
from pathlib import Path
from urllib.parse import unquote

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

# The GenBank LOCUS line puts the identifier and the sequence length inside a
# 28 character field, so the identifier budget shrinks as contigs get longer.
LOCUS_FIELD_WIDTH = 28

# GFF3 feature types worth carrying over, mapped to their INSDC feature keys.
# `gene` and `region` are dropped: they duplicate the CDS/RNA spans and neither
# consumer looks at them.
GFF_TYPE_TO_INSDC = {
    "cds": "CDS",
    "trna": "tRNA",
    "rrna": "rRNA",
    "tmrna": "tmRNA",
    "ncrna": "ncRNA",
    "misc_rna": "misc_RNA",
    "regulatory_region": "regulatory",
    "crispr": "repeat_region",
    "repeat_region": "repeat_region",
    "oric": "rep_origin",
    "oriv": "rep_origin",
    "orit": "rep_origin",
    "gap": "assembly_gap",
    "assembly_gap": "assembly_gap",
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", required=True, type=Path, help="Nucleotide assembly FASTA (authoritative contig IDs).")
    parser.add_argument("--gbk", type=Path, help="Existing GenBank annotation to harmonise.")
    parser.add_argument("--gff", type=Path, help="GFF3 annotation, used when no GenBank file is given.")
    parser.add_argument("--faa", type=Path, help="Protein FASTA supplying CDS translations.")
    parser.add_argument("--output", required=True, type=Path, help="Output GenBank file.")
    parser.add_argument("--contig-map", required=True, type=Path, help="Output contig identifier mapping table.")
    parser.add_argument("--organism", default="", help="Organism name written to the GenBank SOURCE field.")
    parser.add_argument("--transl-table", type=int, default=11, help="Translation table used when a CDS has no translation.")
    args = parser.parse_args(argv)

    if not args.gbk and not args.gff:
        parser.error("one of --gbk or --gff is required")
    return args


def smart_open(path: Path):
    """Open plain or gzip-compressed text, transparently."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_fasta(path: Path) -> "OrderedDict[str, SeqRecord]":
    with smart_open(path) as handle:
        records = OrderedDict()
        for record in SeqIO.parse(handle, "fasta"):
            records[record.id] = record
    if not records:
        sys.exit(f"ERROR: no sequences found in {path}")
    return records


def read_protein_fasta(path: Path | None) -> dict[str, str]:
    """Map protein ID (first header token, as written by Prokka and Bakta) to sequence."""
    if not path:
        return {}
    with smart_open(path) as handle:
        # Internal stops would make the GenBank /translation unparseable for
        # strict consumers, so they are masked the same way as in translate_cds.
        return {
            record.id: str(record.seq).rstrip("*").replace("*", "X")
            for record in SeqIO.parse(handle, "fasta")
        }


def strip_version(identifier: str) -> str:
    return re.sub(r"\.\d+$", "", identifier)


def safe_locus(identifier: str, seq_length: int) -> str:
    """Shorten an identifier so it fits the GenBank LOCUS field.

    Truncated names keep a hash of the original so that distinct contigs cannot
    collapse onto the same LOCUS.
    """
    sanitised = re.sub(r"[^A-Za-z0-9_.\-]", "_", identifier)
    budget = LOCUS_FIELD_WIDTH - 1 - len(str(seq_length))
    if len(sanitised) <= max(budget, 16):
        return sanitised
    digest = hashlib.md5(identifier.encode()).hexdigest()[:6]
    return f"{sanitised[: max(budget, 16) - 7]}_{digest}"


def translate_cds(sequence: Seq, transl_table: int) -> str:
    """Translate a CDS, tolerating the partial genes common in draft assemblies."""
    trimmed = sequence[: len(sequence) - (len(sequence) % 3)]
    if not trimmed:
        return ""
    return str(trimmed.translate(table=transl_table)).rstrip("*").replace("*", "X")


def parse_gff_attributes(field: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for part in field.split(";"):
        if "=" not in part:
            continue
        key, _, value = part.partition("=")
        attributes[key.strip()] = unquote(value.strip())
    return attributes


def features_from_gff(path: Path) -> dict[str, list[SeqFeature]]:
    """Read GFF3 features, stopping at the ``##FASTA`` section Prokka and Bakta append."""
    features: dict[str, list[SeqFeature]] = {}
    with smart_open(path) as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            seqid, _source, gff_type, start, end, _score, strand, phase, attribute_field = fields[:9]
            insdc_type = GFF_TYPE_TO_INSDC.get(gff_type.lower())
            if not insdc_type:
                continue

            attributes = parse_gff_attributes(attribute_field)
            # GFF3 is 1-based inclusive, Biopython locations are 0-based half-open.
            location = SimpleLocation(
                int(start) - 1,
                int(end),
                strand={"+": 1, "-": -1}.get(strand, 0),
            )

            qualifiers: dict[str, list[str]] = {}
            locus_tag = attributes.get("locus_tag") or attributes.get("ID")
            if locus_tag:
                qualifiers["locus_tag"] = [locus_tag]
            gene_name = attributes.get("gene") or attributes.get("Name")
            if gene_name:
                qualifiers["gene"] = [gene_name]
            qualifiers["product"] = [attributes.get("product", "hypothetical protein")]
            if attributes.get("Dbxref"):
                qualifiers["db_xref"] = attributes["Dbxref"].split(",")
            if insdc_type == "CDS":
                qualifiers["codon_start"] = [str(int(phase) + 1) if phase.isdigit() else "1"]

            features.setdefault(seqid, []).append(SeqFeature(location=location, type=insdc_type, qualifiers=qualifiers))
    return features


def resolve_contig(record: SeqRecord, fasta_records: "OrderedDict[str, SeqRecord]") -> str | None:
    """Match a GenBank record to a FASTA contig via its accession or LOCUS name."""
    versionless = {strip_version(key): key for key in fasta_records}
    for candidate in (record.id, record.name, strip_version(record.id), strip_version(record.name)):
        if not candidate:
            continue
        if candidate in fasta_records:
            return candidate
        if candidate in versionless:
            return versionless[candidate]
    return None


def features_from_genbank(
    path: Path, fasta_records: "OrderedDict[str, SeqRecord]"
) -> tuple[dict[str, list[SeqFeature]], list[str]]:
    features: dict[str, list[SeqFeature]] = {}
    unmatched: list[str] = []
    with smart_open(path) as handle:
        for record in SeqIO.parse(handle, "genbank"):
            contig = resolve_contig(record, fasta_records)
            if contig is None:
                unmatched.append(record.id)
                continue
            kept = [feature for feature in record.features if feature.type != "source"]
            features.setdefault(contig, []).extend(kept)
    return features, unmatched


def finalise_cds(
    feature: SeqFeature,
    record_sequence: Seq,
    protein_sequences: dict[str, str],
    transl_table: int,
) -> None:
    """Ensure a CDS carries the translation IslandPath-DIMOB needs to build its .faa."""
    if feature.qualifiers.get("translation"):
        return

    locus_tag = (feature.qualifiers.get("locus_tag") or feature.qualifiers.get("protein_id") or [""])[0]
    translation = protein_sequences.get(locus_tag)
    if not translation:
        translation = translate_cds(feature.extract(record_sequence), transl_table)
    if translation:
        feature.qualifiers["translation"] = [translation]
    feature.qualifiers.setdefault("transl_table", [str(transl_table)])


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    fasta_records = read_fasta(args.fasta)
    protein_sequences = read_protein_fasta(args.faa)

    if args.gbk:
        features_by_contig, unmatched = features_from_genbank(args.gbk, fasta_records)
        for record_id in unmatched:
            print(
                f"WARNING: GenBank record '{record_id}' has no matching contig in {args.fasta.name} and was dropped",
                file=sys.stderr,
            )
    else:
        features_by_contig = features_from_gff(args.gff)
        for seqid in set(features_by_contig) - set(fasta_records):
            print(
                f"WARNING: GFF sequence '{seqid}' has no matching contig in {args.fasta.name} and was dropped",
                file=sys.stderr,
            )

    output_records: list[SeqRecord] = []
    mapping_rows: list[dict[str, object]] = []

    for contig_id, fasta_record in fasta_records.items():
        sequence = fasta_record.seq
        features = features_by_contig.get(contig_id, [])

        for feature in features:
            if feature.type == "CDS":
                finalise_cds(feature, sequence, protein_sequences, args.transl_table)

        locus = safe_locus(contig_id, len(sequence))
        record = SeqRecord(
            sequence,
            id=contig_id,
            name=locus,
            description=fasta_record.description[len(fasta_record.id) :].strip() or contig_id,
            features=features,
        )
        record.annotations["molecule_type"] = "DNA"
        record.annotations["topology"] = "linear"
        record.annotations["organism"] = args.organism or "Bacteria"
        record.annotations["source"] = args.organism or "Bacteria"
        output_records.append(record)

        mapping_rows.append(
            {
                "contig_id": contig_id,
                "locus_id": locus,
                "accession_id": contig_id,
                "length": len(sequence),
                "n_cds": sum(1 for feature in features if feature.type == "CDS"),
                "n_features": len(features),
            }
        )

    total_cds = sum(int(row["n_cds"]) for row in mapping_rows)
    if total_cds == 0:
        sys.exit(
            "ERROR: no CDS features could be assigned to any contig. IslandPath-DIMOB and PhiSpy "
            "both require annotated CDS features with a matching nucleotide sequence."
        )

    with open(args.output, "w") as handle:
        SeqIO.write(output_records, handle, "genbank")

    with open(args.contig_map, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(mapping_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(mapping_rows)

    print(f"Wrote {len(output_records)} contigs and {total_cds} CDS features to {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
