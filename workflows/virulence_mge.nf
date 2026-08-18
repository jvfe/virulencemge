/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { PREPARE_GENOMES        } from '../subworkflows/local/prepare_genomes'
include { VIRULENCE_SCREENING    } from '../subworkflows/local/virulence_screening'
include { MGE_PREDICTION         } from '../subworkflows/local/mge_prediction'
include { ELEMENT_INTEGRATION    } from '../subworkflows/local/element_integration'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VIRULENCE_MGE {

    take:
    ch_samplesheet // channel: [ val(meta), path(fasta), path(gff), path(gbk), path(faa) ]
    outdir         //  string: output directory

    main:

    def ch_versions = channel.empty()

    //
    // SUBWORKFLOW: Harmonise annotation files into one GenBank file per genome
    //
    PREPARE_GENOMES ( ch_samplesheet )

    //
    // SUBWORKFLOW: Virulence factor mining with ABricate against VFDB
    //
    VIRULENCE_SCREENING (
        PREPARE_GENOMES.out.fasta,
        params.abricate_datadir ? file(params.abricate_datadir, checkIfExists: true) : [],
    )

    //
    // SUBWORKFLOW: Prophage and genomic island prediction
    //
    // Nextflow's strict syntax parser hands every command line value over as a
    // string, so `--skip_phispy false` would otherwise be truthy.
    MGE_PREDICTION (
        PREPARE_GENOMES.out.gbk,
        params.skip_phispy.toString().toBoolean(),
        params.skip_islandpath.toString().toBoolean(),
    )

    //
    // SUBWORKFLOW: Coordinate integration, cross-referencing and reporting
    //
    ELEMENT_INTEGRATION (
        PREPARE_GENOMES.out.contig_map,
        VIRULENCE_SCREENING.out.report,
        MGE_PREDICTION.out.prophages,
        MGE_PREDICTION.out.islands,
        VIRULENCE_SCREENING.out.summary,
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'virulencemge_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    emit:
    annotation = ELEMENT_INTEGRATION.out.annotation // channel: [ val(meta), path(tsv) ]
    report     = ELEMENT_INTEGRATION.out.report     // channel: path(html)
    versions   = ch_collated_versions               // channel: path(versions.yml)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
