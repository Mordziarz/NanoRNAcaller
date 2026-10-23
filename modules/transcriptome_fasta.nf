process TRANSCRIPTOME_FASTA {
    tag "Transcriptome FASTA"
    label 'process_medium'

    publishDir "${params.outdir}/transcriptome", mode: 'copy', overwrite: true

    input:
    path reference_genome
    path gtf

    output:
    path "transcriptome.fasta", emit: fasta

    script:
    """
    set -euo pipefail

    echo "Generating transcriptome sequence using gffread..."
    gffread -w transcriptome.fasta -g ${reference_genome} ${gtf}

    echo "Transcriptome successfully generated."
    """
}
