process GET_DORADO {
    tag "Get latest Dorado"
    label 'process_low'

    publishDir "${params.outdir}/dorado",
        mode: 'symlink',
        overwrite: true

    output:
    path "dorado_bin", emit: bin_dir

    script:
    """
    set -euo pipefail

    echo "Checking latest Dorado release on GitHub..."

    REDIRECT_URL=\$(curl -sSI https://github.com/nanoporetech/dorado/releases/latest | grep -i "^location:" | tail -n 1 | tr -d '\\r' | awk '{print \$2}')

    if [ -z "\$REDIRECT_URL" ]; then
        echo "ERROR: Failed to determine redirect to the latest release." >&2
        exit 1
    fi

    RAW_VERSION=\$(basename "\$REDIRECT_URL")
    VERSION=\$(echo "\$RAW_VERSION" | sed 's/^v//')

    echo "Detected latest Dorado version: \$VERSION"

    DOWNLOAD_URL="https://cdn.oxfordnanoportal.com/software/analysis/dorado-\$VERSION-linux-x64.tar.gz"

    echo "Downloading from official CDN: \$DOWNLOAD_URL"

    curl -fSL --retry 3 --connect-timeout 30 \
        "\$DOWNLOAD_URL" -o dorado_archive

    mkdir extracted
    tar -xzf dorado_archive -C extracted

    EXTRACTED_DIR=\$(find extracted \
        -mindepth 1 \
        -maxdepth 1 \
        -type d \
        -name 'dorado-*-linux-x64*' \
        -print -quit)

    if [ -z "\$EXTRACTED_DIR" ]; then
        echo "ERROR: Dorado directory not found after extraction." >&2
        find extracted -maxdepth 3 -print
        exit 1
    fi

    rm -rf dorado_bin
    mv "\$EXTRACTED_DIR" dorado_bin

    echo "Verifying installation:"
    ./dorado_bin/bin/dorado --version

    rm -rf extracted dorado_archive
    echo "Dorado \$VERSION has been successfully downloaded and configured."
    """
}
