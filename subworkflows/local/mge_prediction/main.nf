//
// Predict the mobile genetic elements that virulence factors may reside in
//
// Both tools read the harmonised GenBank file produced by PREPARE_GENOMES:
// PhiSpy calls integrated prophages, IslandPath-DIMOB calls genomic islands from
// dinucleotide bias next to mobility genes.
//

include { PHISPY     } from '../../../modules/nf-core/phispy'
include { ISLANDPATH } from '../../../modules/nf-core/islandpath'

workflow MGE_PREDICTION {

    take:
    ch_gbk          // channel: [ val(meta), path(gbk) ]
    skip_phispy     // boolean: skip prophage prediction
    skip_islandpath // boolean: skip genomic island prediction

    main:

    ch_prophages            = channel.empty()
    ch_prophage_annotations = channel.empty()
    ch_prophage_gff         = channel.empty()
    ch_islands              = channel.empty()

    if (!skip_phispy) {
        PHISPY ( ch_gbk )
        ch_prophages            = PHISPY.out.coordinates
        ch_prophage_annotations = PHISPY.out.information
        ch_prophage_gff         = PHISPY.out.prophage_gff
    }

    if (!skip_islandpath) {
        ISLANDPATH ( ch_gbk )
        ch_islands = ISLANDPATH.out.gff
    }

    emit:
    prophages            = ch_prophages            // channel: [ val(meta), path(tsv)  ]
    prophage_annotations = ch_prophage_annotations  // channel: [ val(meta), path(tsv)  ]
    prophage_gff         = ch_prophage_gff          // channel: [ val(meta), path(gff3) ]
    islands              = ch_islands               // channel: [ val(meta), path(gff)  ]
}
