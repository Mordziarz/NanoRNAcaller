process APALORD {
    tag "APALORD_Analysis"
    label 'process_high'

    publishDir "${params.outdir}/apalord", mode: 'copy'

    input:
    path isoquant_dir
    path isoquant_gtf
    path samples_csv
    val is_drna

    output:
    path "apalord_out/*", emit: results

    script:
    def drna_flag = is_drna ? "--direct_rna TRUE" : "--direct_rna FALSE"
    """
    Rscript ${projectDir}/bin/run_apalord.R \
        --isoquant_dir ${isoquant_dir} \
        --gtf ${isoquant_gtf} \
        --samples_csv ${samples_csv} \
        --cores ${task.cpus} \
        ${drna_flag} \
        --outdir apalord_out
    """
}
