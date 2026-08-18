process INTEGRATE_ELEMENTS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_pandas:1d6d425241d4db42'
        : 'community.wave.seqera.io/library/python_pandas:29a3c74bdab36f58'}"

    input:
    tuple val(meta), path(abricate), path(phispy), path(islandpath), path(contig_map)

    output:
    tuple val(meta), path("${prefix}.virulence_annotation.tsv"), emit: annotation
    tuple val(meta), path("${prefix}.mge_regions.tsv")         , emit: regions
    tuple val(meta), path("${prefix}.stats.tsv")               , emit: stats
    tuple val(meta), path("${prefix}.mge_regions.bed")         , emit: regions_bed
    tuple val(meta), path("${prefix}.mge_regions.gff3")        , emit: regions_gff
    tuple val(meta), path("${prefix}.virulence_factors.bed")   , emit: virulence_bed
    tuple val("${task.process}"), val('pandas'), eval('python -c "import pandas; print(pandas.__version__)"'), topic: versions, emit: versions_pandas

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def abricate_arg = abricate ? "--abricate ${abricate}" : ''
    def phispy_arg = phispy ? "--phispy ${phispy}" : ''
    def islandpath_arg = islandpath ? "--islandpath ${islandpath}" : ''
    def contig_map_arg = contig_map ? "--contig-map ${contig_map}" : ''
    """
    integrate_elements.py \\
        --sample ${meta.id} \\
        --prefix ${prefix} \\
        ${abricate_arg} \\
        ${phispy_arg} \\
        ${islandpath_arg} \\
        ${contig_map_arg} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.virulence_annotation.tsv
    touch ${prefix}.mge_regions.tsv
    touch ${prefix}.stats.tsv
    touch ${prefix}.mge_regions.bed
    touch ${prefix}.mge_regions.gff3
    touch ${prefix}.virulence_factors.bed
    """
}
