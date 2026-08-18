<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-virulencemge_logo_dark.png">
    <img alt="nf-core/virulencemge" src="docs/images/nf-core-virulencemge_logo_light.png">
  </picture>
</h1>

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/nf-core/virulencemge)
[![GitHub Actions CI Status](https://github.com/nf-core/virulencemge/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nf-core/virulencemge/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/virulencemge/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/virulencemge/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/virulencemge/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.1.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.1.0)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/virulencemge)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23virulencemge-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/virulencemge)[![Follow on Bluesky](https://img.shields.io/badge/bluesky-%40nf__core-1185fe?labelColor=000000&logo=bluesky)](https://bsky.app/profile/nf-co.re)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/virulencemge** takes bacterial genomes that have already been annotated with
[Prokka](https://github.com/tseemann/prokka) or [Bakta](https://github.com/oschwengers/bakta) and works out
which of their virulence factors sit inside mobile genetic elements. It screens the assemblies against
[VFDB](http://www.mgc.ac.cn/VFs/) with [ABricate](https://github.com/tseemann/abricate), predicts integrated
prophages with [PhiSpy](https://github.com/linsalrob/PhiSpy) and genomic (pathogenicity) islands with
[IslandPath-DIMOB](https://github.com/brinkmanlab/islandpath), then intersects the three sets of coordinates.
The output is a per-genome table stating, for every virulence gene, whether it lies in a pathogenicity island,
a prophage, both, or the core genome, alongside BED/GFF3 tracks and a standalone HTML dashboard.

1. Harmonise the annotation into a single GenBank file per genome, rebuilding it from GFF3 + FASTA when no
   GenBank file was supplied (`prepare_genbank.py`)
2. Mine virulence factors ([`ABricate`](https://github.com/tseemann/abricate), VFDB) and aggregate them into a
   presence/absence matrix (`abricate --summary`)
3. Identify prophages ([`PhiSpy`](https://github.com/linsalrob/PhiSpy))
4. Identify genomic / pathogenicity islands ([`IslandPath-DIMOB`](https://github.com/brinkmanlab/islandpath))
5. Cross-reference the coordinates of all three and classify each virulence factor's genomic context
   (`integrate_elements.py`)
6. Render summary matrices and a self-contained HTML report (`summarise_virulence_mge.py`)

### Why the annotation is rebuilt

PhiSpy and IslandPath-DIMOB both read GenBank, but they report contig identifiers differently: PhiSpy uses the
`ACCESSION` line while IslandPath-DIMOB uses the `LOCUS` line, and either can differ from the FASTA header that
ABricate reports. Intersecting their coordinates naively silently produces zero overlaps. The pipeline therefore
rewrites the GenBank file so all three identifiers derive from the assembly FASTA, and emits a mapping table that
the integration step uses to reconcile whatever the tools report.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet describing your pre-annotated genomes:

`samplesheet.csv`:

```csv
sample,fasta,gff,gbk,faa
MT1012,/data/MT1012.fna,/data/MT1012.gff,/data/MT1012.gbk,/data/MT1012.faa
MT1013,/data/MT1013.fna,/data/MT1013.gff3,/data/MT1013.gbff,
MT1014,/data/MT1014.fna.gz,/data/MT1014.gff3.gz,,
```

`sample` and `fasta` are required. `gff`, `gbk` and `faa` are optional columns that may be left empty, but each
row needs at least one of `gbk` or `gff`: when `gbk` is missing the pipeline builds one from `gff` + `fasta`, and
uses `faa` for the CDS translations if it is available. Any file may be gzip-compressed.

Now, you can run the pipeline using:

```bash
nextflow run nf-core/virulencemge \
   -profile <docker/singularity/conda> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/virulencemge/usage) and the [parameter documentation](https://nf-co.re/virulencemge/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/virulencemge/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/virulencemge/output).

## Credits

nf-core/virulencemge was originally written by jvfe.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#virulencemge` channel](https://nfcore.slack.com/channels/virulencemge) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

If you use nf-core/virulencemge for your analysis, please cite the pipeline release you used together with
the tools listed in [`CITATIONS.md`](CITATIONS.md).

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
