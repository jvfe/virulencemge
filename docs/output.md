# nf-core/virulencemge: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Prepared genomes](#prepared-genomes) - Harmonised GenBank annotation and contig identifier map
- [ABricate](#abricate) - Virulence factor detection against VFDB
- [PhiSpy](#phispy) - Prophage identification
- [IslandPath-DIMOB](#islandpath-dimob) - Genomic / pathogenicity island identification
- [Integrated elements](#integrated-elements) - Per-genome coordinate cross-referencing and visualisation tracks
- [Summary report](#summary-report) - Multi-genome matrices and standalone HTML dashboard
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Prepared genomes

<details markdown="1">
<summary>Output files</summary>

- `prepared_genomes/`
  - `*_prepared.gbk`: GenBank annotation used as input to PhiSpy and IslandPath-DIMOB. Either your own GenBank file rewritten so that its `LOCUS` and `ACCESSION` identifiers come from the assembly FASTA, or a new file built from the GFF3 and FASTA. Nucleotide sequence and CDS `/translation` qualifiers are always present, because IslandPath-DIMOB fails without them.
  - `*_prepared.contig_map.tsv`: Mapping from FASTA contig identifier to the `LOCUS` and `ACCESSION` identifiers written into the GenBank file, plus contig length and CDS count.

</details>

This step exists because the three analysis tools disagree about how to name contigs. PhiSpy reports the GenBank
`ACCESSION`, IslandPath-DIMOB reports the `LOCUS` name (with the version suffix stripped) and ABricate reports the
FASTA header. The contig map is what lets the integration step line the three coordinate systems up; if any
identifier still cannot be reconciled it is counted in `n_unresolved_contig_ids` and flagged in the HTML report.

### ABricate

<details markdown="1">
<summary>Output files</summary>

- `abricate/`
  - `<sample>.txt`: Tabular [ABricate](https://github.com/tseemann/abricate) report of virulence factor hits for one genome.
  - `abricate_summary.txt`: Presence/absence matrix across all genomes, produced by `abricate --summary`. Columns are genes, cells hold the percentage coverage of the hit.

</details>

The database is selected with `--abricate_db` (default `vfdb`) and hit stringency with `--abricate_minid` and
`--abricate_mincov`.

### PhiSpy

<details markdown="1">
<summary>Output files</summary>

- `phispy/`
  - `<sample>.tsv`: Prophage coordinates (`prophage_coordinates.tsv`), the file used for cross-referencing. Columns are prophage number, contig, start, stop and, where detected, the _attL_/_attR_ site coordinates and sequences.
  - `<sample>.gbk`: Input GenBank annotated with the predicted prophage regions.
  - `<sample>.log`: PhiSpy run log.
  - `<sample>_prophage.tsv`, `<sample>_prophage.tbl`, `<sample>_prophage.gff3`: Simplified prophage region tables and GFF3 track.
  - `<sample>_prophage_information.tsv`: Per-gene table with the phage-likeness score of every gene; the tenth column is `0` for genes considered bacterial and otherwise holds the prophage number. This is the file to inspect when assessing whether a prophage call is trustworthy.
  - `<sample>_phage.gbk`, `<sample>_bacteria.gbk`: Prophage and prophage-free portions of the genome.

</details>

Which of these are produced is controlled by `--phispy_output_choice` (default `512`, meaning all of them).

### IslandPath-DIMOB

<details markdown="1">
<summary>Output files</summary>

- `islandpath/`
  - `<sample>.gff`: Predicted genomic islands in GFF3. Islands are called from regions of atypical dinucleotide bias that also contain a mobility gene (integrase, transposase and similar).
  - `logs/<sample>_Dimob.log`: IslandPath-DIMOB run log.

</details>

IslandPath-DIMOB takes no tuning options; the minimum island size is fixed at 8 kb by the tool's own
configuration file.

### Integrated elements

<details markdown="1">
<summary>Output files</summary>

- `integrated_elements/`
  - `<sample>.virulence_annotation.tsv`: One row per virulence factor hit. `mge_context` is the headline column and holds `PAI`, `Prophage`, `PAI+Prophage` or `Core genome`. `pai_ids`/`prophage_ids` name the containing elements, and `*_overlap_bp`/`*_overlap_fraction` quantify the overlap so borderline calls can be reviewed.
  - `<sample>.mge_regions.tsv`: One row per predicted element with its coordinates, the virulence genes it carries, and `cross_element_ids`/`cross_element_overlap_bp` describing overlap with elements of the other type. Genomic islands frequently contain prophages, and this is where that shows up.
  - `<sample>.stats.tsv`: Single-row count summary for the genome, including `frac_vf_in_mge` and `n_unresolved_contig_ids`.
  - `<sample>.mge_regions.bed`, `<sample>.mge_regions.gff3`: Prophage and island tracks for genome browsers such as IGV, Artemis or JBrowse.
  - `<sample>.virulence_factors.bed`: Virulence factor track; each feature is named `<gene>|<context>`.

</details>

A virulence factor counts as contained in an element when at least `--mge_min_overlap` of the gene (default `0.5`)
falls inside it. Lower the threshold to call any overlap, or set it to `1` to require full containment.

### Summary report

<details markdown="1">
<summary>Output files</summary>

- `summary_report/`
  - `virulencemge_report.html`: Standalone dashboard. Plotly is inlined, so the file needs no network access. Contains summary cards, quality warnings, virulence factor burden and proportion plots, element counts and size distributions, the co-localisation frequency plot, a gene-by-genome localisation heatmap and the underlying tables.
  - `virulencemge_sample_summary.tsv`: Per-genome counts, concatenated across all samples.
  - `virulencemge_vf_localisation_matrix.tsv`: Gene by genome matrix whose cells state the genomic context (`PAI`, `Prophage`, `PAI+Prophage`, `Core genome` or `Absent`). Where a gene is present in several copies in one genome, the most informative context wins, so an MGE-borne copy is never masked by a core-genome copy.
  - `virulencemge_vf_presence_matrix.tsv`: The same matrix reduced to 0/1 presence.
  - `virulencemge_virulence_annotation.tsv`, `virulencemge_mge_regions.tsv`: All per-genome tables concatenated.
  - `virulencemge_abricate_summary.tsv`: Copy of the `abricate --summary` matrix.

</details>

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
