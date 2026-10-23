#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Rsamtools)
  library(dplyr)
  library(stats)
  library(tidyr)
  library(rlang)
  library(ggplot2)
  library(ComplexHeatmap)
  library(ggtree)
  library(circlize)
  library(car)
  library(dunn.test)
  library(nortest)
  library(rstatix)
  library(data.table)
  library(rtracklayer)
  library(TAILcaller)
})

args <- commandArgs(trailingOnly = TRUE)
samples_csv <- args[1]
gtf_file_path <- args[2]

cat("[INFO] Reading samples table...\n")
df_samples <- read.csv(samples_csv, stringsAsFactors = FALSE)

# TAILcaller expects columns: bam_path, sample_name, group
samples_table <- data.frame(
  bam_path = df_samples$bam_path,
  sample_name = df_samples$sample_name,
  group = df_samples$group,
  stringsAsFactors = FALSE
)

cat("[INFO] Running TAILcaller::get_polyA()...\n")
get_polyA_out <- TAILcaller::get_polyA(samples_table = samples_table)

cat("[INFO] Filtering valid poly(A) lengths (> 0)...\n")
get_polyA_out <- get_polyA_out[get_polyA_out$polyA_length > 0, ]
fwrite(get_polyA_out, "polyA_raw_lengths.tsv", sep = "\t")

cat("[INFO] Importing GTF file and assigning gene IDs...\n")
gtf_raw <- import.gff(gtf_file_path)
gtf_df <- as.data.frame(gtf_raw)

get_gene_id_out <- TAILcaller::get_gene_id(
  polyA_table = get_polyA_out,
  gtf_file = gtf_df,
  transcript_column_gtf = "transcript_id",
  gene_column_gtf = "gene_id"
)
fwrite(get_gene_id_out, "polyA_with_gene_id.tsv", sep = "\t")

unique_groups <- unique(samples_table$group)
if (length(unique_groups) >= 2) {
  ctrl_g <- unique_groups[1]
  treat_g <- unique_groups[2]
  cat(sprintf("[INFO] Running statistics comparing control: '%s' vs treated: '%s'\n", ctrl_g, treat_g))

  tryCatch({
    cat("[INFO] Calculating transcript-level statistics (calculate_statistics)...\n")
    stat_tx_1 <- TAILcaller::calculate_statistics(
      polyA_table = get_gene_id_out,
      grouping_factor = "group",
      which_level = "transcript_id",
      control_group = ctrl_g,
      treated_group = treat_g,
      padj_method = "fdr"
    )
    fwrite(stat_tx_1, "tail_stats_transcript_level_1.tsv", sep = "\t")
  }, error = function(e) cat("[WARN] Transcript stats 1 error:", e$message, "\n"))

  tryCatch({
    cat("[INFO] Calculating transcript-level statistics (calculate_statistics_n2)...\n")
    stat_tx_2 <- TAILcaller::calculate_statistics_n2(
      polyA_table = get_gene_id_out,
      grouping_factor = "group",
      which_level = "transcript_id",
      control_group = ctrl_g,
      treated_group = treat_g,
      padj_method = "fdr"
    )
    fwrite(stat_tx_2, "tail_stats_transcript_level_2.tsv", sep = "\t")
  }, error = function(e) cat("[WARN] Transcript stats 2 error:", e$message, "\n"))

  tryCatch({
    cat("[INFO] Calculating gene-level statistics (calculate_statistics)...\n")
    stat_gene_1 <- TAILcaller::calculate_statistics(
      polyA_table = get_gene_id_out,
      grouping_factor = "group",
      which_level = "gene_id",
      control_group = ctrl_g,
      treated_group = treat_g,
      padj_method = "fdr"
    )
    fwrite(stat_gene_1, "tail_stats_gene_level_1.tsv", sep = "\t")
  }, error = function(e) cat("[WARN] Gene stats 1 error:", e$message, "\n"))

  tryCatch({
    cat("[INFO] Calculating gene-level statistics (calculate_statistics_n2)...\n")
    stat_gene_2 <- TAILcaller::calculate_statistics_n2(
      polyA_table = get_gene_id_out,
      grouping_factor = "group",
      which_level = "gene_id",
      control_group = ctrl_g,
      treated_group = treat_g,
      padj_method = "fdr"
    )
    fwrite(stat_gene_2, "tail_stats_gene_level_2.tsv", sep = "\t")
  }, error = function(e) cat("[WARN] Gene stats 2 error:", e$message, "\n"))
}

cat("[INFO] TAILcaller analysis script completed successfully.\n")
