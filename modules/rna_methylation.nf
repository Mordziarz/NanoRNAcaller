process RNA_METHYLATION {
    tag "RNA_Methylation_${sample} (${assay})"
    label 'process_high'

    publishDir "${params.outdir}/rna_methylation/${sample}", mode: 'copy', pattern: "*.{bam,csi,bed,log}"

    input:
    tuple val(sample), val(group), val(assay), path(aligned_bam), path(aligned_index)
    path ref_genome

    output:
    tuple val(sample), val(group), val(assay), path(aligned_bam), path(aligned_index), path("${sample}_nuc.bed"), emit: methylation_results
    path "${sample}_pileup.log", emit: log

    script:
    """
    set -eo pipefail

    echo "[INFO] Running modkit pileup on pre-sorted BAM for ${sample}..."
    modkit pileup \\
        --ref ${ref_genome} \\
        --log-filepath ${sample}_pileup.log \\
        ${aligned_bam} \\
        ${sample}_nuc.bed

    echo "[INFO] RNA methylation analysis completed for ${sample}."
    """
}
