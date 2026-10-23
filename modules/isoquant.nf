process ISOQUANT {
    tag "IsoQuant_Analysis"
    label 'process_high'

    publishDir "${params.outdir}/isoquant", mode: 'copy'

    input:
    val samples
    path bams
    path bais
    val assays
    path reference_genome
    path gtf

    output:
    path "isoquant_out", emit: isoquant_dir
    path "isoquant_out/**/*", emit: results

    script:
    def sample_labels = samples.join(' ')
    def bam_files     = bams.join(' ')
    
    def is_drna = assays.any { it.toLowerCase() == 'drs' || it.toLowerCase().contains('drna') }
    def stranded_arg = is_drna ? "--stranded forward" : ""

    """
    isoquant \
        --reference ${reference_genome} \
        --genedb ${gtf} \
        --large_output read_assignments corrected_bed \
        -d nanopore \
        ${stranded_arg} \
        --labels ${sample_labels} \
        --bam ${bam_files} \
        -o isoquant_out
        
    find isoquant_out -name "*.gz" -exec gunzip {} +
    """
}
