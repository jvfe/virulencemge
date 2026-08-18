//
// Screen assemblies for virulence factors and aggregate them into one matrix
//

include { ABRICATE_RUN     } from '../../../modules/nf-core/abricate/run'
include { ABRICATE_SUMMARY } from '../../../modules/nf-core/abricate/summary'

workflow VIRULENCE_SCREENING {

    take:
    ch_fasta     // channel: [ val(meta), path(fasta) ]
    ch_databases //   value: path to a custom ABricate database directory, or []

    main:

    ABRICATE_RUN ( ch_fasta, ch_databases )

    //
    // `abricate --summary` turns the per-genome reports into a single
    // presence/absence matrix, so it has to see every report at once.
    //
    ABRICATE_SUMMARY (
        ABRICATE_RUN.out.report
            .map { _meta, report -> report }
            .collect()
            .map { reports -> [ [ id: 'abricate_summary' ], reports ] }
    )

    emit:
    report  = ABRICATE_RUN.out.report      // channel: [ val(meta), path(txt) ]
    summary = ABRICATE_SUMMARY.out.report  // channel: [ val(meta), path(txt) ]
}
