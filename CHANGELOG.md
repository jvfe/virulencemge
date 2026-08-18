# nf-core/virulencemge: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of nf-core/virulencemge, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Samplesheet input taking pre-annotated bacterial genomes (`sample,fasta,gff,gbk,faa`), where the annotation may
  be supplied as GenBank or as GFF3.
- `PREPARE_GENBANK`, building a GenBank file with contig identifiers, nucleotide sequence and CDS translations that
  satisfy the requirements of the downstream tools.
- Virulence factor screening with ABricate against VFDB, aggregated across genomes with `abricate --summary`.
- Prophage prediction with PhiSpy and genomic island prediction with IslandPath-DIMOB, each skippable.
- `INTEGRATE_ELEMENTS`, cross-referencing virulence factor coordinates against predicted prophages and genomic
  islands, and emitting per-genome tables plus BED and GFF3 tracks.
- `SUMMARY_REPORT`, producing multi-genome matrices and a self-contained HTML dashboard.

### `Fixed`

### `Dependencies`

| Dependency   | Version |
| ------------ | ------- |
| `abricate`   | 1.0.1   |
| `phispy`     | 4.2.21  |
| `islandpath` | 1.0.6   |
| `biopython`  | 1.84    |
| `pandas`     | 2.2.3   |
| `plotly`     | 5.24.1  |
| `jinja2`     | 3.1.4   |

### `Deprecated`
