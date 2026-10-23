#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(edgeR)
    library(ggplot2)
    library(reshape2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    stop("Usage: run_transcript_de_analysis.R <transcript_counts_file> <samples_table_file> [output_dir]", call. = FALSE)
}

counts_file <- args[1]
meta_file   <- args[2]
out_dir     <- if (length(args) >= 3) args[3] else "transcript_de_results"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

counts <- read.delim(counts_file, header = TRUE, row.names = 1, check.names = FALSE)
meta   <- read.csv(meta_file, header = TRUE)

common_samples <- intersect(colnames(counts), meta$sample)
if (length(common_samples) == 0) {
    colnames(counts) <- gsub(".sorted.bam", "", colnames(counts))
    common_samples <- intersect(colnames(counts), meta$sample)
}

counts <- counts[, common_samples, drop = FALSE]
meta   <- meta[match(common_samples, meta$sample), ]

group_counts <- table(meta$group)
min_replicates <- min(group_counts)
cat("Minimum number of replicates per group (Transcripts):", min_replicates, "\n")

cpm_vals <- cpm(counts, log = TRUE)
pca <- prcomp(t(cpm_vals), scale. = TRUE)
pca_df <- as.data.frame(pca$x)
pca_df$sample <- rownames(pca_df)
pca_df$group  <- meta$group

var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group, label = sample)) +
    geom_point(size = 4) +
    geom_text(vjust = -1, size = 3, show.legend = FALSE) +
    theme_minimal() +
    labs(
        title = "PCA - Global Transcript Expression",
        x = paste0("PC1 (", var_explained[1], "%)"),
        y = paste0("PC2 (", var_explained[2], "%)")
    )

ggsave(file.path(out_dir, "transcript_pca_plot.pdf"), p_pca, width = 7, height = 5)
ggsave(file.path(out_dir, "transcript_pca_plot.png"), p_pca, width = 7, height = 5, dpi = 300)

# 3. edgeR analysis (always performed)
y <- DGEList(counts = counts, group = meta$group)
y <- calcNormFactors(y)
design <- model.matrix(~ 0 + meta$group)
colnames(design) <- levels(factor(meta$group))

y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

groups <- levels(factor(meta$group))
if (length(groups) == 2) {
    contrast_formula <- makeContrasts(contrasts = paste0(groups[1], "-", groups[2]), levels = design)
    qlf <- glmQLFTest(fit, contrast = contrast_formula)
    edger_res <- topTags(qlf, n = Inf)$table
    
    # Dołączenie surowych zliczeń oraz znormalizowanych CPM do tabeli edgeR dla transkryptów
    raw_counts_subset <- counts[rownames(edger_res), , drop = FALSE]
    colnames(raw_counts_subset) <- paste0("raw_", colnames(raw_counts_subset))
    
    cpm_subset <- cpm(y)[rownames(edger_res), , drop = FALSE]
    colnames(cpm_subset) <- paste0("cpm_norm_", colnames(cpm_subset))
    
    edger_final <- cbind(edger_res, raw_counts_subset, cpm_subset)
    write.csv(edger_final, file.path(out_dir, "edger_differential_transcript_expression.csv"))
}

if (min_replicates >= 3) {
    cat("All groups have >= 3 replicates. Running DESeq2 for transcripts...\n")
    suppressPackageStartupMessages(library(DESeq2))
    
    colData <- data.frame(row.names = meta$sample, group = factor(meta$group))
    dds <- DESeqDataSetFromMatrix(
        countData = round(counts),
        colData   = colData,
        design    = ~ group
    )
    dds <- DESeq(dds)
    
    if (length(groups) == 2) {
        res <- results(dds)
        res_df <- as.data.frame(res)
        
        norm_counts_deseq <- counts(dds, normalized = TRUE)[rownames(res_df), , drop = FALSE]
        colnames(norm_counts_deseq) <- paste0("deseq_norm_", colnames(norm_counts_deseq))
        
        raw_counts_deseq <- counts(dds, normalized = FALSE)[rownames(res_df), , drop = FALSE]
        colnames(raw_counts_deseq) <- paste0("raw_", colnames(raw_counts_deseq))
        
        deseq_final <- cbind(res_df, raw_counts_deseq, norm_counts_deseq)
        write.csv(deseq_final, file.path(out_dir, "deseq2_differential_transcript_expression.csv"))
    }
} else {
    cat("WARNING: At least one group has fewer than 3 replicates (min:", min_replicates, "). Skipping DESeq2 for transcripts.\n")
    file.create(file.path(out_dir, "DESeq2_transcript_skipped_due_to_low_replicates.txt"))
}
