process MINIMAP2_ALIGN {
    tag "Align_${sample} (${assay})"
    label 'process_high'

    publishDir "${params.outdir}/alignment", mode: 'copy', pattern: "*.sorted.bam*"

    input:
    tuple val(sample), val(group), val(assay), path(fastq_file)
    path ref_genome

    output:
    tuple val(sample), val(group), val(assay), path("${sample}.sorted.bam"), path("${sample}.sorted.bam.csi"), emit: aligned_bam

    script:
    def minimap2_opts = (assay == "drs" || assay.contains("drna")) ? "-ax splice -uf -k14 -w4 -ub" : "-ax splice -k14 -w4 -ub"

    """
    set -eo pipefail

    echo "[INFO] Building minimap2 index with increased capacity (-I 8g)..."
    minimap2 -t ${task.cpus} -I 8g -d reference.mmi ${ref_genome}

    echo "[INFO] Running alignment and sorting..."
    minimap2 -t ${task.cpus} ${minimap2_opts} reference.mmi ${fastq_file} | \
    samtools view -b -h -F 2308 - | \
    samtools sort -@ ${task.cpus} -o ${sample}.sorted.bam -

    echo "[INFO] Indexing BAM file..."
    samtools index -c ${sample}.sorted.bam
    """
}
