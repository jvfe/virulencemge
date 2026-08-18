//
// Cross-reference virulence factors against prophages and genomic islands
//
// The per-sample step decides whether each VFDB hit sits inside a mobile genetic
// element; the aggregation step turns those calls into multi-sample matrices and
// a standalone HTML dashboard.
//

include { INTEGRATE_ELEMENTS } from '../../../modules/local/integrate_elements'
include { SUMMARY_REPORT     } from '../../../modules/local/summary_report'

workflow ELEMENT_INTEGRATION {

    take:
    ch_contig_map       // channel: [ val(meta), path(tsv) ]
    ch_abricate         // channel: [ val(meta), path(txt) ]
    ch_prophages        // channel: [ val(meta), path(tsv) ]
    ch_islands          // channel: [ val(meta), path(gff) ]
    ch_abricate_summary // channel: [ val(meta), path(txt) ]

    main:

    //
    // The contig map exists for every sample, so it drives the join; any tool
    // that was skipped or produced nothing contributes a null that the module
    // turns into an omitted argument.
    //
    ch_integration_input = ch_contig_map
        .join( ch_abricate,  remainder: true )
        .join( ch_prophages, remainder: true )
        .join( ch_islands,   remainder: true )
        .map { meta, contig_map, abricate, prophages, islands ->
            [ meta, abricate ?: [], prophages ?: [], islands ?: [], contig_map ]
        }

    INTEGRATE_ELEMENTS ( ch_integration_input )

    SUMMARY_REPORT (
        INTEGRATE_ELEMENTS.out.stats.map { _meta, stats -> stats }.collect(),
        INTEGRATE_ELEMENTS.out.annotation.map { _meta, annotation -> annotation }.collect(),
        INTEGRATE_ELEMENTS.out.regions.map { _meta, regions -> regions }.collect(),
        ch_abricate_summary.map { _meta, summary -> summary }.ifEmpty( [] )
    )

    emit:
    annotation     = INTEGRATE_ELEMENTS.out.annotation     // channel: [ val(meta), path(tsv)  ]
    regions        = INTEGRATE_ELEMENTS.out.regions        // channel: [ val(meta), path(tsv)  ]
    stats          = INTEGRATE_ELEMENTS.out.stats          // channel: [ val(meta), path(tsv)  ]
    regions_bed    = INTEGRATE_ELEMENTS.out.regions_bed    // channel: [ val(meta), path(bed)  ]
    regions_gff    = INTEGRATE_ELEMENTS.out.regions_gff    // channel: [ val(meta), path(gff3) ]
    virulence_bed  = INTEGRATE_ELEMENTS.out.virulence_bed  // channel: [ val(meta), path(bed)  ]
    report         = SUMMARY_REPORT.out.report             // channel: path(html)
    sample_summary = SUMMARY_REPORT.out.sample_summary     // channel: path(tsv)
}
