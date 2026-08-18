# nf-core/virulencemge: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/virulencemge/usage](https://nf-co.re/virulencemge/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

nf-core/virulencemge annotates virulence factors, prophages and genomic (pathogenicity) islands in bacterial
genomes that have **already been annotated**, typically with Prokka or Bakta, and then reports which virulence
factors sit inside mobile genetic elements. It does not perform assembly or gene calling.

## Samplesheet input

You will need to create a samplesheet describing your annotated genomes and pass it with `--input`. It is a
comma-separated file with a header row and the five columns below.

```bash
--input '[path to samplesheet file]'
```

```csv title="samplesheet.csv"
sample,fasta,gff,gbk,faa
MT1012,/data/MT1012.fna,/data/MT1012.gff,/data/MT1012.gbk,/data/MT1012.faa
MT1013,/data/MT1013.fna,/data/MT1013.gff3,/data/MT1013.gbff,
MT1014,/data/MT1014.fna.gz,/data/MT1014.gff3.gz,,
```

| Column   | Description                                                                                                                                                                                                        |
| -------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `sample` | Custom genome identifier. Must be unique and cannot contain spaces; it is used to name every output file.                                                                                                          |
| `fasta`  | **Required.** Nucleotide assembly FASTA (`.fasta`, `.fas`, `.fa`, `.fna`, `.fsa`, optionally `.gz`). This is what ABricate is screened against, and its headers define the contig identifiers used in all outputs. |
| `gff`    | Optional. Annotation GFF3 (`.gff`, `.gff3`, optionally `.gz`), as written by Prokka or Bakta. Used to build a GenBank file when `gbk` is absent.                                                                   |
| `gbk`    | Optional. Annotation GenBank file (`.gbk`, `.gb`, `.gbf`, `.gbff`, optionally `.gz`).                                                                                                                              |
| `faa`    | Optional. Protein FASTA (`.faa`, `.fa`, `.fasta`, optionally `.gz`). Supplies CDS translations when building a GenBank file from GFF3.                                                                             |

### Optional columns

Only `sample` and `fasta` are mandatory, but **every row needs at least one of `gbk` or `gff`**, because PhiSpy and
IslandPath-DIMOB can only read GenBank. Rows are validated individually, so a run can freely mix genomes that have
a GenBank file with genomes that only have a GFF3.

- **`gbk` given** (recommended): it is used directly, after being rewritten so its contig identifiers match the FASTA.
- **`gbk` empty, `gff` given**: a GenBank file is built from `gff` + `fasta`. CDS translations are taken from `faa`
  when present and otherwise generated with translation table 11. Prokka and Bakta both append the assembly to
  their GFF3 under a `##FASTA` header; that section is ignored in favour of the `fasta` column.
- **`faa` empty**: harmless, translations are computed as described above.

Leaving a column empty means an empty field in the CSV (`,,`), not the string `NA`. An
[example samplesheet](../assets/samplesheet.csv) is provided with the pipeline.

### Annotation requirements

The GenBank file handed to PhiSpy and IslandPath-DIMOB must contain CDS features **and** the nucleotide sequence;
IslandPath-DIMOB aborts otherwise. The pipeline guarantees this by injecting the sequence from the `fasta` column
and filling in any missing `/translation` qualifiers, so annotations that would fail if passed to the tools
directly will usually work here. What it cannot fix is a genome with no CDS features at all, which fails with an
explicit error.

## Selecting a database

Virulence factors are screened with ABricate against VFDB by default. Any database bundled with the ABricate
container can be selected instead:

```bash
--abricate_db card
```

To use a locally built or more recent copy of VFDB, point `--abricate_datadir` at your ABricate database directory.
Hit stringency is controlled with `--abricate_minid` and `--abricate_mincov`.

## Tuning the co-localisation call

A virulence factor is reported as residing in an element when at least `--mge_min_overlap` of the gene (default
`0.5`) falls inside the predicted prophage or genomic island. Use a small value to call any overlap at all, or `1`
to require the gene to be fully contained:

```bash
--mge_min_overlap 0.99
```

The per-genome `*.virulence_annotation.tsv` always reports the raw overlap in base pairs and as a fraction, so the
threshold can be revisited without re-running the tools.

## Skipping steps

`--skip_phispy` and `--skip_islandpath` turn off prophage and genomic island prediction respectively. The
integration and reporting steps still run, and virulence factors are simply reported against whichever elements
remain.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/virulencemge --input ./samplesheet.csv --outdir ./results  -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/virulencemge -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/virulencemge
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/virulencemge releases page](https://github.com/nf-core/virulencemge/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
