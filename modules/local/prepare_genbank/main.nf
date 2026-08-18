process PREPARE_GENBANK {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/python_biopython:b26cebed064c3984'
        : 'community.wave.seqera.io/library/python_biopython:fa00db0dce19aedd'}"

    input:
    tuple val(meta), path(fasta), path(gff), path(gbk), path(faa)

    output:
    tuple val(meta), path("${prefix}.gbk")           , emit: gbk
    tuple val(meta), path("${prefix}.contig_map.tsv"), emit: contig_map
    tuple val("${task.process}"), val('biopython'), eval('python -c "import Bio; print(Bio.__version__)"'), topic: versions, emit: versions_biopython

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}_prepared"
    // An existing GenBank file is preferred; GFF3 is only used to rebuild one from scratch
    def annotation = gbk ? "--gbk ${gbk}" : "--gff ${gff}"
    def faa_arg = faa ? "--faa ${faa}" : ''

    if ("${gbk}" == "${prefix}.gbk") {
        error "Input and output GenBank names are the same for '${meta.id}', set prefix in the module configuration to disambiguate!"
    }

    """
    prepare_genbank.py \\
        --fasta ${fasta} \\
        ${annotation} \\
        ${faa_arg} \\
        --output ${prefix}.gbk \\
        --contig-map ${prefix}.contig_map.tsv \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_prepared"
    """
    touch ${prefix}.gbk
    echo -e "contig_id\\tlocus_id\\taccession_id\\tlength\\tn_cds\\tn_features" > ${prefix}.contig_map.tsv
    """
}
