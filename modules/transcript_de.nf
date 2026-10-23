process DIFFERENTIAL_TRANSCRIPT_EXPRESSION {
    tag "Transcript DE (${transcript_counts.baseName})"
    label 'process_medium'

    publishDir "${params.outdir}/transcript_differential_expression", mode: 'copy'

    input:
    path transcript_counts
    path samples_csv
    path r_script

    output:
    path "transcript_de_results/*", emit: results

    script:
    """
    Rscript ${r_script} \
        ${transcript_counts} \
        ${samples_csv} \
        transcript_de_results
    """
}
