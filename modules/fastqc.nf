process FASTQC {
    tag "$meta.id"
    label 'process_medium'
    container 'mordziarz/transcriptomelongpipeline:latest'

    env _JAVA_OPTIONS: "-Xmx16G"

    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"),  emit: zip

    script:
    """
    fastqc -t ${task.cpus} ${reads}
    """
}
