process DIFFERENTIAL_EXPRESSION {
    tag "Diff_Expression_Analysis"
    label 'process_medium'

    publishDir "${params.outdir}/differential_expression", mode: 'copy'

    input:
    path gene_counts
    path samples_table
    path script_file

    output:
    path "de_results/*", emit: results

    script:
    """
    Rscript ${script_file} ${gene_counts} ${samples_table} de_results
    """
}
