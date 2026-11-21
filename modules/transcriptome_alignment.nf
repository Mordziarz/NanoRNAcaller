process TRANSCRIPTOME_ALIGNMENT {
    tag "Transcriptome_Align_${sample} (${assay})"
    label 'process_high'

    publishDir "${params.outdir}/transcriptome_alignment/${sample}", mode: 'copy', pattern: "*.{bam,csi}"
    
    input:
    tuple val(sample), val(group), val(assay), path(bam_file)
    path transcriptome_fasta
    path dorado_bin

    output:
    tuple val(sample),
          val(group),
          val(assay),
          path("${sample}_aligned.sorted.bam"),
          path("${sample}_aligned.sorted.bam.csi"),
          emit: aligned_bam

    script:
    def mm2_options = (assay == "drs" || assay.contains("drna")) ? "-x splice" : "-x map-ont"

    """
    set -eo pipefail

    echo "[INFO] Running dorado transcriptome aligner for sample ${sample} (assay: ${assay})..."

    ${dorado_bin}/bin/dorado aligner \\
        --mm2-opts "${mm2_options}" \\
        ${transcriptome_fasta} \\
        ${bam_file} \\
        > ${sample}_raw.sam

    echo "[INFO] Filtering with samtools (-F 2308) and sorting for ${sample}..."

    samtools view -b -h -F 2308 ${sample}_raw.sam | \\
        samtools sort -@ ${task.cpus} \\
        -o ${sample}_aligned.sorted.bam -

    echo "[INFO] Checking sorted BAM file..."

    samtools quickcheck -v ${sample}_aligned.sorted.bam

    echo "[INFO] Indexing sorted BAM file using CSI..."

    samtools index -c -@ ${task.cpus} \\
        ${sample}_aligned.sorted.bam

    test -s ${sample}_aligned.sorted.bam.csi

    rm -f ${sample}_raw.sam

    echo "[INFO] Dorado transcriptome alignment and filtering completed for ${sample}."
    """
}
