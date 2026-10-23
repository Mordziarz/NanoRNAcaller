process BAM_TO_FASTQ {
    tag "BAM2FASTQ_${sample} (${assay})"
    label 'process_medium'

    publishDir "${params.outdir}/fastq", mode: 'copy', pattern: "*.fastq.gz"

    input:
    tuple val(sample), val(group), val(assay), path(bam_file)

    output:
    tuple val(sample), val(group), val(assay), path("${sample}.fastq.gz"), emit: fastq

    script:
    """
    samtools fastq -T '*' ${bam_file} | gzip > ${sample}.fastq.gz
    """
}
