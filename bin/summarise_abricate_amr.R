# summarise_abricate_amr.R
# Summarise & visualise AMR calls from Abricate reports
# Usage: setwd("path/with/abricate_reports"); source("summarise_abricate_amr.R")

# ---- User settings ----
min_coverage <- 90   # %COVERAGE threshold
min_identity <- 90   # %IDENTITY threshold
out_prefix   <- "amr_summary"

# ---- Libraries ----
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(janitor)
  library(ggplot2)
  library(pheatmap)
  library(forcats)
})

# ---- Load all abricate reports ----
files <- list.files(pattern = "_abricate_report\\.txt$", full.names = TRUE)
if (length(files) == 0) stop("No *_abricate_report.txt files found in working directory.")

read_abricate <- function(f) {
  # Explicit column names per Abricate report format
  coln <- c("FILE","SEQUENCE","START","END","STRAND","GENE","ACCESSION",
            "COVERAGE","GAPS","PCT_COVERAGE","PCT_IDENTITY",
            "DATABASE","DB_ACCESSION","PRODUCT","RESISTANCE")
  
  df <- readr::read_tsv(
    f,
    comment = "#",                # skip the '#...' header lines
    col_names = coln,             # prevent first data row becoming header
    col_types = readr::cols(.default = readr::col_character())
  )
  df$source_file <- basename(f)
  df
}

# ---- Load all reports ----
files <- list.files(pattern = "_abricate_report\\.txt$", full.names = TRUE)
if (length(files) == 0) stop("No *_abricate_report.txt files found.")
raw <- dplyr::bind_rows(lapply(files, read_abricate))

# ---- Clean & standardise columns ----
df <- raw %>%
  janitor::clean_names() %>%
  mutate(
    # Parse numeric from strings like "100.00" (no % sign in Abricate, but safe if present)
    percent_coverage = as.numeric(stringr::str_replace(pct_coverage, "%", "")),
    percent_identity = as.numeric(stringr::str_replace(pct_identity, "%", "")),
    gene         = dplyr::coalesce(gene, ""),
    resistance   = dplyr::coalesce(resistance, ""),
    file_fasta   = dplyr::coalesce(file, source_file),
    # Sample ID from FILE or fallback to filename
    sample = dplyr::case_when(
      stringr::str_detect(file_fasta, "INF\\d+") ~ stringr::str_extract(file_fasta, "INF\\d+"),
      stringr::str_detect(source_file, "INF\\d+") ~ stringr::str_extract(source_file, "INF\\d+"),
      TRUE ~ tools::file_path_sans_ext(source_file)
    )
  )


# ---- Filter on QC thresholds ----
df_filt <- df %>%
  filter(percent_coverage >= min_coverage,
         percent_identity >= min_identity)

# ---- Basic summaries ----
# 1) Genes per sample
genes_per_sample <- df_filt %>%
  count(sample, gene, sort = TRUE)

# 2) Resistance classes per sample (split semicolon-separated values)
res_classes <- df_filt %>%
  mutate(class = if_else(resistance == "", NA_character_, resistance)) %>%
  separate_longer_delim(class, delim = ";") %>%
  mutate(class = str_replace_all(class, "_", " ")) %>%
  filter(!is.na(class)) %>%
  count(sample, class, sort = TRUE)

# 3) Presence/absence matrix (genes x samples)
presence <- genes_per_sample %>%
  mutate(value = 1L) %>%
  select(sample, gene, value) %>%
  distinct() %>%
  pivot_wider(names_from = sample, values_from = value, values_fill = 0) %>%
  arrange(gene)

# ---- Write tables ----
write_csv(df_filt, paste0(out_prefix, "_filtered_calls.csv"))
write_csv(genes_per_sample, paste0(out_prefix, "_genes_per_sample.csv"))
write_csv(res_classes, paste0(out_prefix, "_resistance_classes.csv"))
write_csv(presence, paste0(out_prefix, "_presence_absence_matrix.csv"))

# ---- Plots ----
dir.create("plots", showWarnings = FALSE)

# A) Bar: number of AMR genes per sample
genes_count <- genes_per_sample %>%
  group_by(sample) %>%
  summarise(n_genes = n_distinct(gene), .groups = "drop")

p1 <- ggplot(genes_count, aes(x = fct_reorder(sample, n_genes), y = n_genes)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "AMR genes detected per sample",
    x = "Sample",
    y = "Number of genes (filtered)"
  ) +
  theme_minimal(base_size = 12)

ggsave(filename = file.path("plots", paste0(out_prefix, "_genes_per_sample_bar.png")),
       plot = p1, width = 6, height = 4, dpi = 300)

# B) Heatmap: gene presence/absence
mat <- presence %>%
  column_to_rownames("gene") %>%
  as.matrix()

# Keep informative ordering (cluster by default)
png(file.path("plots", paste0(out_prefix, "_presence_absence_heatmap.png")),
    width = 1100, height = 900, res = 150)
pheatmap(mat,
         main = "AMR gene presence/absence",
         border_color = NA,
         legend = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE)
dev.off()

# C) Stacked bar: resistance classes per sample
# Keep top 10 classes visible; group the rest into 'Other'
top_classes <- res_classes %>%
  group_by(class) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  slice_max(n, n = 10) %>%
  pull(class)

res_classes_plot <- res_classes %>%
  mutate(class_grp = if_else(class %in% top_classes, class, "Other")) %>%
  group_by(sample, class_grp) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(frac = n / sum(n))

p2 <- ggplot(res_classes_plot, aes(x = sample, y = frac, fill = class_grp)) +
  geom_col() +
  labs(
    title = "AMR resistance classes per sample",
    x = "Sample",
    y = "Proportion of calls",
    fill = "Class"
  ) +
  theme_minimal(base_size = 12)

ggsave(filename = file.path("plots", paste0(out_prefix, "_resistance_classes_stacked.png")),
       plot = p2, width = 7, height = 4.5, dpi = 300)

# ---- Console message ----
message("Done.\n",
        "- Tables written: ",
        paste(c(paste0(out_prefix, "_filtered_calls.csv"),
                paste0(out_prefix, "_genes_per_sample.csv"),
                paste0(out_prefix, "_resistance_classes.csv"),
                paste0(out_prefix, "_presence_absence_matrix.csv")), collapse = ", "),
        "\n- Plots in ./plots/: ",
        paste(c(paste0(out_prefix, "_genes_per_sample_bar.png"),
                paste0(out_prefix, "_presence_absence_heatmap.png"),
                paste0(out_prefix, "_resistance_classes_stacked.png")), collapse = ", "),
        "\nThresholds: coverage ≥ ", min_coverage, "%, identity ≥ ", min_identity, "%")