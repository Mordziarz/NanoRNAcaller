process GENOME_ALIGNMENT {
    tag "Genome_Align_${sample} (${assay})"
    label 'process_high'

    input:
    tuple val(sample), val(group), val(assay), path(bam_file)
    path ref_genome
    path dorado_bin

    output:
    tuple val(sample),
          val(group),
          val(assay),
          path("${sample}_genome_aligned.sorted.bam"),
          path("${sample}_genome_aligned.sorted.bam.csi"),
          emit: aligned_bam

    script:
    def mm2_options = (assay == "drs" || assay.contains("drna")) ? "-x splice" : "-x map-ont"

    """
    set -eo pipefail

    echo "[INFO] Running dorado genome aligner for sample ${sample} (assay: ${assay})..."

    ${dorado_bin}/bin/dorado aligner \\
        --mm2-opts "${mm2_options}" \\
        ${ref_genome} \\
        ${bam_file} \\
        > ${sample}_raw.sam

    echo "[INFO] Filtering with samtools (-F 2308) and sorting for ${sample}..."

    samtools view -b -h -F 2308 ${sample}_raw.sam | \\
        samtools sort -@ ${task.cpus} \\
        -o ${sample}_genome_aligned.sorted.bam -

    echo "[INFO] Checking sorted BAM file..."

    samtools quickcheck -v ${sample}_genome_aligned.sorted.bam

    echo "[INFO] Indexing sorted BAM file using CSI..."

    samtools index -c -@ ${task.cpus} \\
        ${sample}_genome_aligned.sorted.bam

    test -s ${sample}_genome_aligned.sorted.bam.csi

    rm -f ${sample}_raw.sam

    echo "[INFO] Dorado genome alignment and filtering completed for ${sample}."
    """
}
