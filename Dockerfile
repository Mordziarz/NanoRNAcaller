FROM nvidia/cuda:12.8.1-base-ubuntu22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    wget \
    git \
    bzip2 \
    jq \
    tar \
    build-essential \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    software-properties-common \
    dirmngr \
    gpg \
    samtools \
    minimap2 \
    gffread \
    fastqc \
    && add-apt-repository ppa:deadsnakes/ppa \
    && apt-get update && apt-get install -y --no-install-recommends \
    python3.12 \
    python3.12-dev \
    python3.12-venv \
    && curl -sS https://bootstrap.pypa.io/get-pip.py | python3.12 \
    && python3.12 -m pip install --no-cache-dir setuptools wheel \
    && rm -rf /var/lib/apt/lists/*

RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.12 1

RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O miniforge.sh && \
    bash miniforge.sh -b -p /opt/conda && \
    rm miniforge.sh

ENV PATH="/opt/conda/bin:$PATH"

RUN wget https://github.com/nanoporetech/modkit/releases/download/v0.6.4/modkit_v0.6.4_u16_x86_64.tar.gz && \
    tar -xzf modkit_v0.6.4_u16_x86_64.tar.gz && \
    find . -name "modkit" -type f -exec cp {} /usr/local/bin/ \; && \
    rm -rf modkit_v0.6.4_u16_x86_64.tar.gz

RUN git clone https://github.com/comprna/SUPPA.git /opt/suppa
RUN sed -i 's/from statsmodels.sandbox.stats.multicomp import multipletests/from statsmodels.stats.multitest import multipletests/g' /opt/suppa/lib/diff_tools.py

RUN mamba install -y -c conda-forge -c bioconda \
    python=3.12 \
    isoquant \
    statsmodels \
    r-base=4.5.2 \
    r-recommended \
    r-remotes \
    r-data.table \
    r-optparse \
    r-ggplot2 \
    r-reshape2 \
    r-pbmcapply \
    r-stringr \
    r-dplyr \
    r-tidyr \
    r-purrr \
    r-car \
    r-dunn.test \
    r-nortest \
    r-rstatix \
    r-circlize \
    bioconductor-edger \
    bioconductor-deseq2 \
    bioconductor-dexseq \
    bioconductor-genomicranges \
    bioconductor-rtracklayer \
    bioconductor-iranges \
    bioconductor-biostrings \
    bioconductor-genomicalignments \
    bioconductor-rsamtools \
    bioconductor-genomicfeatures \
    bioconductor-bambu \
    bioconductor-complexheatmap \
    bioconductor-ggtree \
    && mamba clean -a -y

RUN git clone https://github.com/markandtwin/APALORD.git /opt/apalord && \
    R CMD INSTALL /opt/apalord

RUN git clone https://github.com/Mordziarz/TAILcaller.git /opt/tailcaller && \
    R CMD INSTALL /opt/tailcaller

RUN pip3 install --no-cache-dir torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu124
RUN pip3 install --no-cache-dir scipy pysam numpy pandas scikit-learn progressbar whatshap

WORKDIR /data
