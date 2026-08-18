//
// Normalise the annotation files supplied in the samplesheet
//
// PhiSpy and IslandPath-DIMOB only read GenBank, and IslandPath-DIMOB
// additionally needs the nucleotide sequence and CDS translations to be present
// in that file. This subworkflow guarantees a single uncompressed GenBank file
// per sample whose contig identifiers agree with the assembly FASTA that
// ABricate is screened against.
//

include { GUNZIP           } from '../../../modules/nf-core/gunzip'
include { PREPARE_GENBANK  } from '../../../modules/local/prepare_genbank'

workflow PREPARE_GENOMES {

    take:
    ch_samplesheet // channel: [ val(meta), path(fasta), path(gff), path(gbk), path(faa) ]

    main:

    //
    // ABricate is the only tool reading the FASTA directly, so that is the one
    // file we decompress here; the GenBank builder reads gzipped input itself.
    //
    ch_input = ch_samplesheet
        .branch { _meta, fasta, _gff, _gbk, _faa ->
            compressed:   fasta.toString().endsWith('.gz')
            uncompressed: true
        }

    GUNZIP (
        ch_input.compressed.map { meta, fasta, _gff, _gbk, _faa -> [ meta, fasta ] }
    )

    ch_genomes = GUNZIP.out.gunzip
        .join( ch_input.compressed.map { meta, _fasta, gff, gbk, faa -> [ meta, gff, gbk, faa ] } )
        .mix( ch_input.uncompressed )

    PREPARE_GENBANK ( ch_genomes )

    emit:
    fasta      = ch_genomes.map { meta, fasta, _gff, _gbk, _faa -> [ meta, fasta ] } // channel: [ val(meta), path(fasta) ]
    gbk        = PREPARE_GENBANK.out.gbk                                             // channel: [ val(meta), path(gbk) ]
    contig_map = PREPARE_GENBANK.out.contig_map                                      // channel: [ val(meta), path(tsv) ]
}
