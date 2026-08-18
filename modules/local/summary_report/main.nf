process SUMMARY_REPORT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_pandas_plotly_jinja2:f7c62a80c7461368'
        : 'community.wave.seqera.io/library/python_pandas_plotly_jinja2:2f8e0791b9d1050d'}"

    input:
    path stats       , stageAs: 'stats/*'
    path annotations , stageAs: 'annotations/*'
    path regions     , stageAs: 'regions/*'
    path abricate_summary

    output:
    path "${prefix}_report.html"        , emit: report
    path "${prefix}_sample_summary.tsv" , emit: sample_summary
    path "${prefix}_*.tsv"              , emit: tables
    tuple val("${task.process}"), val('plotly'), eval('python -c "import plotly; print(plotly.__version__)"'), topic: versions, emit: versions_plotly

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: 'virulencemge'
    def regions_arg = regions ? "--elements regions/*" : ''
    def abricate_arg = abricate_summary ? "--abricate-summary ${abricate_summary}" : ''
    """
    summarise_virulence_mge.py \\
        --stats stats/* \\
        --annotations annotations/* \\
        ${regions_arg} \\
        ${abricate_arg} \\
        --prefix ${prefix} \\
        --pipeline-version '${workflow.manifest.version}' \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: 'virulencemge'
    """
    touch ${prefix}_report.html
    touch ${prefix}_sample_summary.tsv
    touch ${prefix}_vf_localisation_matrix.tsv
    """
}
