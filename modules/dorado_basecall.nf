process DORADO_BASECALL {
    tag "Basecalling_${sample} (${assay})"
    label 'process_high'

    input:
    tuple val(sample), val(group), val(assay), path(pod5_file)
    path dorado_bin

    output:
    tuple val(sample), val(group), val(assay), path("${sample}.bam"), emit: bam

    script:
    def extra_args = (assay == "drs" || assay.contains("drna")) ? 
        "--modified-bases inosine_m6A_2OmeA pseU_2OmeU m5C_2OmeC 2OmeG --emit-moves --estimate-poly-a" : 
        "--emit-moves --trim all --estimate-poly-a"

    """
    set -eo pipefail

    DORADO_BIN="\$(pwd)/dorado_bin/bin/dorado"

    echo "[INFO] Running basecalling (auto-detecting sup model) for ${sample} (assay: ${assay})..."
    
    "\${DORADO_BIN}" basecaller sup ${extra_args} ${pod5_file} > ${sample}.bam
    
    echo "[INFO] Basecalling completed successfully for ${sample}."
    """
}
