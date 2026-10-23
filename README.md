# NanoRNAcaller: A Nextflow Pipeline for Analyzing Direct RNA and cDNA Reads

[![Nextflow](https://img.shields.io/badge/Nextflow-%E2%89%A521.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/Docker-containerized-blue.svg)](https://www.docker.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## Overview

**NanoRNAcaller** is a scalable, containerized Nextflow pipeline designed for end-to-end processing of Oxford Nanopore Technologies (ONT) direct RNA (DRS) and complementary DNA (cDNA) long-read transcriptomics data. 

Unlike traditional analytical workflows that rely heavily on pretrained models, species-specific gene models, or predefined untranslated region (UTR) boundaries, NanoRNAcaller is completely **reference-agnostic**. It operates efficiently on any organism provided with a reference genome assembly and basic annotation file, making it ideal for both model and non-model organisms.

---

## Key Features

* **Dual-Protocol Support:** Hardcoded processing logic optimized for both direct RNA sequencing (DRS) and cDNA datasets.
* **Automated Basecalling:** Fetches the latest stable release of Dorado automatically and executes basecalling using the super accuracy (`sup`) model.
* **Epitranscriptomic Modification Detection:** Identifies and quantifies methylation marks directly from raw DRS signal files using Modkit.
* **Transcript Discovery & Quantification:** Performs gene and transcript quantification along with novel transcript identification via IsoQuant.
* **Differential Expression Analysis:** Automated statistical workflows for DGE and DTE using edgeR and DESeq2 based on group sample sizes.
* **Alternative Splicing & Differential transcript usage:** Evaluates alternative splicing (AS) events and differential transcript usage (DTU) using SUPPA2.
* **Alternative Polyadenylation:** Identifies and profiles APA sites across conditions using APAlord.
* **Poly(A) Tail Profiling:** Measures poly(A) tail dynamics and differential tail length using TAILcaller.
* **Zero-Dependency Overhead:** Fully containerized via Docker for single-command execution and maximum reproducibility.

---

## Workflow Architecture

NanoRNAcaller seamlessly integrates the following toolsuite into a unified Nextflow workflow:

1. **Basecalling & Signal Alignment:** Dorado
2. **Quality Control:** FastQC v0.11.9
3. **Read Alignment & Filtering:** Minimap2 v2.31, Samtools v1.24
4. **Transcript Annotation & Extraction:** IsoQuant v4.0.0, GffRead v0.12.7
5. **Modification Calling (DRS only):** Modkit, Samtools v1.24
6. **Differential Expression:** edgeR v4.10.7, DESeq2 v1.50.2, ggplot2 v4.0.3
7. **Alternative Splicing & DTU:** SUPPA2 v2.3
8. **Alternative Polyadenylation:** APAlord v1.0.0
9. **Poly(A) Tail Dynamics:** TAILcaller v1.0.0

---

## Inputs and Requirements

### 1. Prerequisites
* **Nextflow** ($\ge$ 21.04.0)
* **Docker**

### 2. Reference Files
* Reference Genome Assembly (`.fasta`)
* Genome Annotation (`.gtf`)

### 3. Sample Metadata Sheet
A sample table (`samples.csv`) must be provided with the following column structure:

| sample | group | assay | pod5 | bam |
| :--- | :--- | :--- | :--- | :--- |
| Sample1 | Control | drs | `/path/to/sample1.pod5` | `/path/to/sample1.bam` |
| Sample2 | Control | drs | `/path/to/sample2.pod5` | |
| Sample3 | Treated | cdna | `/path/to/sample3.pod5` | `/path/to/sample3.bam` |

* `sample`: Unique identifier for each sample.
* `group`: Biological condition/group required for downstream statistical testing.
* `assay`: Sequencing protocol, either `drs` or `cdna`.
* `pod5`: Path to raw POD5 signal files.
* `bam`: Path to basecalled BAM files (optional). If omitted or empty, NanoRNAcaller automatically initiates basecalling from POD5 files using Dorado.

---

## Parameter Flags Reference

NanoRNAcaller utilizes pre-configured, optimized parameters for internal tool execution:

* **Dorado Basecalling:**
  * **DRS:** `--modified-bases inosine_m6A_2OmeA pseU_2OmeU m5C_2OmeC 2OmeG --emit-moves --estimate-poly-a`
  * **cDNA:** `--emit-moves --trim all --estimate-poly-a`
* **Minimap2 Alignment:**
  * **DRS:** `-ax splice -uf -k14 -w4 -ub`
  * **cDNA:** `-ax splice -k14 -w4 -ub`
* **Samtools Filtering:** `-F 2308`
* **Differential Expression Group Thresholds:**
  * **edgeR:** Triggered when $\ge 1$ sample per group is present.
  * **DESeq2:** Triggered when $\ge 3$ samples per group are present.

---

## Quick Start

Nextflow can fetch and run **NanoRNAcaller** directly from GitHub without manual cloning or downloading.

Run the pipeline directly via remote repository:

```bash
nextflow run Mordziarz/NanoRNAcaller \
  --samples_table samples.csv \
  --reference_genome reference_genome.fasta \
  --gtf annotation.gtf
```

```
results/
├── alignment/
├── apalord/
├── differential_expression/
├── dorado/
├── fastq/
├── fastqc/
├── isoquant/
├── rna_methylation/
├── suppa/
├── tailcaller/
├── transcript_differential_expression/
├── transcriptome/
└── transcriptome_alignment/
```

## Citation & Acknowledgements

If you use **NanoRNAcaller** in your research, please cite the manuscript:

> paper in preparation

### Community Acknowledgements
NanoRNAcaller relies on an ecosystem of open-source bioinformatics software. Users are strongly encouraged to cite the original publications for the integrated tools:
- **Nextflow:** [Di Tommaso et al., 2017](https://doi.org/10.1038/nbt.3820)
- **Dorado & Modkit:** [Oxford Nanopore Technologies](https://github.com/nanoporetech)
- **FastQC:** [Andrews S.](https://github.com/s-andrews/fastqc)
- **Samtools:** [Li et al., 2009](https://doi.org/10.1093/bioinformatics/btp352)
- **Minimap2:** [Li H., 2018](https://doi.org/10.1093/bioinformatics/bty191)
- **IsoQuant:** [Prjibelski et al., 2023](https://doi.org/10.1038/s41587-022-01565-y)
- **GffRead:** [Pertea & Pertea, 2020](https://doi.org/10.12688/f1000research.23297.1)
- **edgeR:** [Robinson et al., 2010](https://doi.org/10.1093/bioinformatics/btp616)
- **DESeq2:** [Love et al., 2014](https://doi.org/10.1186/s13059-014-0550-8)
- **SUPPA2:** [Trincado et al., 2018](https://doi.org/10.1186/s13059-018-1417-1)
- **APAlord:** [Zhang et al., 2026](https://doi.org/10.1038/s41467-026-76377-8)
- **TAILcaller:** [Maździarz et al., 2025](https://doi.org/10.1093/bioadv/vbaf235)

---

## Funding

This study was financially supported by:
* **The National Science Centre, Kraków, Poland:** Grant No. `2024/53/N/NZ8/02829`
* **Foundation for Polish Science (FNP):** START stipend awarded to Mateusz Maździarz.

### Support & Contact
If you have any questions, encounter any issues, or need assistance with the pipeline, please feel free to reach out:
- **Mateusz Maździarz** — Email: [mateusz.mazdziarz@uwm.edu.pl](mailto:mateusz.mazdziarz@uwm.edu.pl)