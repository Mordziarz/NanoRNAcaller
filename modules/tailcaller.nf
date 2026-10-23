process TAILCALLER {
    tag "TAILcaller Analysis"
    label 'process_medium'

    publishDir "${params.outdir}/tailcaller", mode: 'copy', overwrite: true

    input:
    val samples
    val groups
    path bams
    path gtf
    path r_script

    output:
    path "*.tsv", emit: results
    path "tailcaller_samples_table.csv", emit: samples_csv

    script:
    """
    cat << 'EOF' > tailcaller_samples_table.csv
sample_name,group,bam_path
${samples.withIndex().collect { s, i -> "${s},${groups[i]},${bams[i]}" }.join('\n')}
EOF

    Rscript ${r_script} tailcaller_samples_table.csv ${gtf}
    """
}
