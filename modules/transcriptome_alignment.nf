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

    # Krok 1: Dopasowanie odczytów do transkryptomu
    ${dorado_bin}/bin/dorado aligner \\
        --mm2-opts "${mm2_options}" \\
        ${transcriptome_fasta} \\
        ${bam_file} \\
        > ${sample}_raw.sam

    echo "[INFO] Filtering with samtools (-F 2308) and sorting for ${sample}..."

    # Krok 2: Usunięcie secondary/supplementary alignments i sortowanie
    samtools view -b -h -F 2308 ${sample}_raw.sam | \\
        samtools sort -@ ${task.cpus} \\
        -o ${sample}_aligned.sorted.bam -

    echo "[INFO] Checking sorted BAM file..."

    # Krok 3: Podstawowa kontrola integralności BAM
    samtools quickcheck -v ${sample}_aligned.sorted.bam

    echo "[INFO] Indexing sorted BAM file using CSI..."

    # Krok 4: Indeksowanie w formacie CSI
    samtools index -c -@ ${task.cpus} \\
        ${sample}_aligned.sorted.bam

    # Krok 5: Sprawdzenie, czy indeks CSI istnieje
    test -s ${sample}_aligned.sorted.bam.csi

    # Krok 6: Usunięcie pliku tymczasowego
    rm -f ${sample}_raw.sam

    echo "[INFO] Dorado transcriptome alignment and filtering completed for ${sample}."
    """
}
