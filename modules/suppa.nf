process SUPPA_GENERATE_EVENTS {

    tag "SUPPA_GENERATE_EVENTS_${format}"
    label 'process_low'

    publishDir "${params.outdir}/suppa/events", mode: 'copy'

    input:
    path gtf
    val format
    val events

    output:
    path "*.ioe", optional: true, emit: ioe_files
    path "*.ioi", optional: true, emit: ioi_files

    script:
    def eventArg = format == 'ioe' ? "-e ${events instanceof List ? events.join(' ') : events}" : ""

    """
    python3 /opt/suppa/suppa.py generateEvents \
        -i ${gtf} \
        -o local_events \
        -f ${format} \
        ${eventArg}
    """
}

process SUPPA_PSI {

    tag "SUPPA_PSI_${annotation_file.baseName}"
    label 'process_medium'

    publishDir "${params.outdir}/suppa/psi", mode: 'copy'

    input:
    path annotation_file
    path expression_file
    val mode

    output:
    path "*.psi", emit: psi_file

    script:
    """
    python3 ${projectDir}/bin/suppa_clean_expression.py \
        --expression ${expression_file} \
        --output clean_expression.txt

    prefix=\$(basename "${annotation_file}" .ioe)

    python3 /opt/suppa/suppa.py psiPerEvent \
        --ioe-file ${annotation_file} \
        --expression-file clean_expression.txt \
        -o psi_\${prefix}
    """
}


process SUPPA_PSI_ISOFORM {

    tag "SUPPA_PSI_ISOFORM"
    label 'process_medium'

    publishDir "${params.outdir}/suppa/isoform_psi", mode: 'copy'

    input:
    path gtf
    path expression_file

    output:
    path "*.psi", emit: psi_file

    script:
    """
    python3 ${projectDir}/bin/suppa_clean_expression.py \
        --expression ${expression_file} \
        --output clean_expression.txt

    python3 /opt/suppa/suppa.py psiPerIsoform \
        -g ${gtf} \
        -e clean_expression.txt \
        -o isoform_psi
    """
}


process SUPPA_DIFF_SPLICE {

    tag "SUPPA_DIFF_SPLICE_${kind}_${prefix}"
    label 'process_medium'

    publishDir "${params.outdir}/suppa/diff_${kind}", mode: 'copy'

    input:
    tuple val(prefix), path(annotation_file), path(psi_file)
    path expression_file
    path samples_csv
    val kind

    output:
    path "*.dpsi", optional: true, emit: dpsi_file
    path "*.psivec", optional: true, emit: psivec_file

    script:
    """
    python3 ${projectDir}/bin/suppa_clean_expression.py \
        --expression ${expression_file} \
        --output clean_expression.txt

    python3 ${projectDir}/bin/suppa_split_conditions.py \
        --psi ${psi_file} \
        --expression clean_expression.txt \
        --samples ${samples_csv} \
        --outdir split

    python3 /opt/suppa/suppa.py diffSplice \
        --method empirical \
        --input ${annotation_file} \
        --psi split/group1.psi split/group2.psi \
        --tpm split/group1.tpm split/group2.tpm \
        --lower-bound 0.05 \
        --gene-correction \
        -o diff_${kind}_${prefix}
    """
}
