nextflow.enable.dsl=2

include { GET_DORADO                           } from './modules/get_dorado'
include { DORADO_BASECALL                      } from './modules/dorado_basecall'
include { BAM_TO_FASTQ                         } from './modules/bam_to_fastq'
include { FASTQC                         } from './modules/fastqc'
include { MINIMAP2_ALIGN                       } from './modules/minimap2_align'
include { ISOQUANT                             } from './modules/isoquant'
include { DIFFERENTIAL_EXPRESSION              } from './modules/differential_expression'
include { DIFFERENTIAL_TRANSCRIPT_EXPRESSION } from './modules/transcript_de'
include { TRANSCRIPTOME_FASTA                  } from './modules/transcriptome_fasta'
include { TRANSCRIPTOME_ALIGNMENT              } from './modules/transcriptome_alignment'
include { GENOME_ALIGNMENT                     } from './modules/genome_alignment'
include { RNA_METHYLATION                      } from './modules/rna_methylation'
include { TAILCALLER                           } from './modules/tailcaller'
include {
    SUPPA_GENERATE_EVENTS as SUPPA_GENERATE_IOE
    SUPPA_GENERATE_EVENTS as SUPPA_GENERATE_IOI
    SUPPA_PSI
    SUPPA_PSI_ISOFORM
    SUPPA_DIFF_SPLICE as SUPPA_DIFF_EVENTS
    SUPPA_DIFF_SPLICE as SUPPA_DIFF_ISOFORM
} from './modules/suppa'
include { APALORD } from './modules/apalord'


workflow {

    // ============================================================
    // INPUT FILES
    // ============================================================

    ch_ref_genome = file(
        params.reference_genome,
        checkIfExists: true
    )

    ch_gtf = file(
        params.gtf,
        checkIfExists: true
    )

    ch_samples_csv = file(
        params.samples_table,
        checkIfExists: true
    )


    // ============================================================
    // READ SAMPLES TABLE
    // ============================================================

    ch_input = Channel
        .fromPath(params.samples_table, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->

            def sample = row.sample
            def group  = row.group ?: 'unassigned'
            def assay  = row.assay ?: 'cdna'

            def pod5 = row.pod5?.trim()
                ? file(row.pod5.trim(), checkIfExists: true)
                : null

            def bam = row.bam?.trim()
                ? file(row.bam.trim(), checkIfExists: true)
                : null

            return tuple(sample, group, assay, pod5, bam)
        }


    // ============================================================
    // BASECALLING
    // ============================================================

    ch_needs_basecalling = ch_input
        .filter { row ->
            row[3] != null && row[4] == null
        }
        .map { row ->
            tuple(row[0], row[1], row[2], row[3])
        }


    ch_pod5_and_bam = ch_input
        .filter { row ->
            row[3] != null && row[4] != null
        }
        .map { row ->
            tuple(row[0], row[1], row[2], row[4])
        }


    ch_only_bam = ch_input
        .filter { row ->
            row[3] == null && row[4] != null
        }
        .map { row ->
            tuple(row[0], row[1], row[2], row[4])
        }

    GET_DORADO()
    
    DORADO_BASECALL(ch_needs_basecalling, GET_DORADO.out.bin_dir)


    ch_ready_bams = ch_only_bam
        .mix(ch_pod5_and_bam)
        .mix(DORADO_BASECALL.out.bam)


    // ============================================================
    // 1. GENOME ALIGNMENT (DORADO)
    // ============================================================

    ch_drs_ready_bams = ch_ready_bams
        .filter { sample, group, assay, bam ->
            def a = assay?.toLowerCase() ?: ''
            a == 'drs' || a.contains('drna')
        }

    GENOME_ALIGNMENT(
        ch_drs_ready_bams,
        ch_ref_genome,
        GET_DORADO.out.bin_dir
    )

    RNA_METHYLATION(
        GENOME_ALIGNMENT.out.aligned_bam,
        ch_ref_genome
    )


    // ============================================================
    // BAM TO FASTQ
    // ============================================================

    BAM_TO_FASTQ(ch_ready_bams)

    fastqc_input = BAM_TO_FASTQ.out.fastq.map { id, group, type, reads -> 
        def meta = [ id: id, group: group, type: type ]
        return [ meta, reads ]
    }

    FASTQC( fastqc_input )
    
    // ============================================================
    // ALIGNMENT (GENOME - MINIMAP2)
    // ============================================================

    MINIMAP2_ALIGN(
        BAM_TO_FASTQ.out.fastq,
        ch_ref_genome
    )


    ch_mapped_streams = MINIMAP2_ALIGN.out.aligned_bam
        .map { sample, group, assay, bam, bai ->
            tuple(sample, bam, bai, assay)
        }


    ch_samples = ch_mapped_streams
        .map { it[0] }
        .collect()

    ch_bams = ch_mapped_streams
        .map { it[1] }
        .collect()

    ch_bais = ch_mapped_streams
        .map { it[2] }
        .collect()

    ch_assays = ch_mapped_streams
        .map { it[3] }
        .collect()


    // ============================================================
    // ISOQUANT
    // ============================================================

    ISOQUANT(
        ch_samples,
        ch_bams,
        ch_bais,
        ch_assays,
        ch_ref_genome,
        ch_gtf
    )


    // ============================================================
    // TRANSCRIPTOME FASTA GENERATION (GFFREAD)
    // ============================================================

    ch_isoquant_gtf = ISOQUANT.out.results
        .flatten()
        .filter {
            it.toString().endsWith('extended_annotation.gtf')
        }

    ch_isoquant_gtf_val = ch_isoquant_gtf.first()

    TRANSCRIPTOME_FASTA(
        ch_ref_genome,
        ch_isoquant_gtf_val
    )


    // ============================================================
    // TRANSCRIPTOME ALIGNMENT
    // ============================================================

    TRANSCRIPTOME_ALIGNMENT(
        ch_ready_bams,
        TRANSCRIPTOME_FASTA.out.fasta,
        GET_DORADO.out.bin_dir
    )

    ch_transcriptome_aligned_bams = TRANSCRIPTOME_ALIGNMENT.out.aligned_bam


    // ============================================================
    // DIFFERENTIAL GENE EXPRESSION
    // ============================================================

    ch_gene_counts = ISOQUANT.out.results
        .flatten()
        .filter {
            it.toString().endsWith(
                'discovered_gene_grouped_file_name_counts.tsv'
            ) &&
            !it.toString().contains('.linear.')
        }


    DIFFERENTIAL_EXPRESSION(
        ch_gene_counts,
        ch_samples_csv,
        file(
            "${projectDir}/bin/run_de_analysis.R",
            checkIfExists: true
        )
    )

    // ============================================================
    // DIFFERENTIAL TRANSCRIPT EXPRESSION
    // ============================================================

    ch_transcript_counts = ISOQUANT.out.results
        .flatten()
        .filter {
            it.toString().endsWith(
                'discovered_transcript_grouped_file_name_counts.tsv'
            ) &&
            !it.toString().contains('.linear.')
        }

    DIFFERENTIAL_TRANSCRIPT_EXPRESSION(
        ch_transcript_counts,
        ch_samples_csv,
        file(
            "${projectDir}/bin/run_transcript_de_analysis.R",
            checkIfExists: true
        )
    )


    // ============================================================
    // ISOQUANT OUTPUTS FOR SUPPA2 & TAILCALLER
    // ============================================================

    ch_transcript_tpm = ISOQUANT.out.results
        .flatten()
        .filter {
            it.toString().endsWith(
                'discovered_transcript_grouped_file_name_tpm.tsv'
            ) &&
            !it.toString().contains('.linear.')
        }

    ch_transcript_tpm_val = ch_transcript_tpm.first()


    // ============================================================
    // SUPPA2: GENERATE EVENT DEFINITIONS (.ioe)
    // ============================================================

    SUPPA_GENERATE_IOE(
        ch_isoquant_gtf_val,
        'ioe',
        ['SE', 'SS', 'MX', 'RI', 'FL']
    )


    ch_all_ioe = SUPPA_GENERATE_IOE.out.ioe_files
        .flatten()
        .filter {
            it.toString().endsWith('.ioe')
        }


    // ============================================================
    // SUPPA2: GENERATE TRANSCRIPT DEFINITIONS (.ioi)
    // ============================================================

    SUPPA_GENERATE_IOI(
        ch_isoquant_gtf_val,
        'ioi',
        ''
    )


    ch_all_ioi = SUPPA_GENERATE_IOI.out.ioi_files
        .flatten()
        .filter {
            it.toString().endsWith('.ioi')
        }


    // ============================================================
    // SUPPA2: PSI PER EVENT
    // ============================================================

    ch_psi_events = SUPPA_PSI(
        ch_all_ioe,
        ch_transcript_tpm_val,
        'event'
    )


    // ============================================================
    // SUPPA2: PSI PER ISOFORM
    // ============================================================

    ch_psi_isoforms = SUPPA_PSI_ISOFORM(
        ch_isoquant_gtf_val,
        ch_transcript_tpm_val
    )


    // ============================================================
    // SUPPA2: MATCH IOE WITH EVENT PSI
    // ============================================================

    ch_event_psi_keyed = ch_psi_events.psi_file
        .flatten()
        .map { psi ->
            def prefix = psi.baseName.replaceFirst(/^psi_/, '')
            tuple(prefix, psi)
        }


    ch_event_diff_input = ch_all_ioe
        .map { ioe ->
            tuple(ioe.baseName, ioe)
        }
        .join(ch_event_psi_keyed)
        .map { prefix, ioe, psi ->
            tuple(prefix, ioe, psi)
        }


    // ============================================================
    // SUPPA2: DIFFERENTIAL SPLICING
    // ============================================================

    SUPPA_DIFF_EVENTS(
        ch_event_diff_input,
        ch_transcript_tpm_val,
        ch_samples_csv,
        'events'
    )


    // ============================================================
    // SUPPA2: MATCH IOI WITH ISOFORM PSI
    // ============================================================

    ch_isoform_diff_input = ch_all_ioi
        .map { ioi ->
            tuple(ioi.baseName, ioi)
        }
        .combine(ch_psi_isoforms.psi_file)
        .map { prefix, ioi, psi ->
            tuple(prefix, ioi, psi)
        }


    // ============================================================
    // SUPPA2: DIFFERENTIAL ISOFORM USAGE
    // ============================================================

    SUPPA_DIFF_ISOFORM(
        ch_isoform_diff_input,
        ch_transcript_tpm_val,
        ch_samples_csv,
        'isoform_usage'
    )


    // ============================================================
    // APALORD: ALTERNATIVE POLYADENYLATION
    // ============================================================

    def has_drna = false

    try {

        def csvFile = new File(params.samples_table)
        def lines = csvFile.readLines()

        if (lines.size() > 0) {

            def header = lines[0].split(',', -1)*.trim()
            def assayIdx = header.indexOf('assay')

            if (assayIdx != -1) {

                for (int i = 1; i < lines.size(); i++) {

                    def cols = lines[i].split(',', -1)*.trim()

                    if (cols.size() > assayIdx) {

                        def val = cols[assayIdx].toLowerCase()

                        if (val == 'drs' || val.contains('drna')) {
                            has_drna = true
                            break
                        }
                    }
                }

            } else {

                has_drna = lines.any {
                    it.toLowerCase().contains('drs') ||
                    it.toLowerCase().contains('drna')
                }
            }
        }

    } catch (Exception e) {

        has_drna = false
    }


    APALORD(
        ISOQUANT.out.isoquant_dir,
        ch_isoquant_gtf,
        ch_samples_csv,
        has_drna
    )


    // ============================================================
    // TAILCALLER
    // ============================================================

    ch_tailcaller_input = ch_transcriptome_aligned_bams
        .map { row ->
            tuple(row[0], row[1], row[3]) // row[0] = sample, row[1] = group, row[3] = bam
        }

    TAILCALLER(
        ch_tailcaller_input.map { it[0] }.collect(),
        ch_tailcaller_input.map { it[1] }.collect(), 
        ch_tailcaller_input.map { it[2] }.collect(), 
        ch_isoquant_gtf_val,
        file("${projectDir}/bin/run_tailcaller.R", checkIfExists: true)
    )
}
