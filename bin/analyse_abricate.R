#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(optparse)
})

# ---------- CLI ----------
option_list <- list(
  make_option(c("-i","--input"), type="character", default="abricate_merged_all.tsv",
              help="Path to merged ABRicate TSV [default: %default]"),
  make_option(c("-o","--outdir"), type="character", default="abricate_plots",
              help="Output directory for plots [default: %default]"),
  make_option(c("-n","--topN"), type="integer", default=30,
              help="Top N AMR genes to plot in distribution [default: %default]"),
  make_option(c("--gene_column"), type="character", default="PRODUCT",
              help="Which field to treat as the AMR gene label in ResFinder (PRODUCT or GENE) [default: %default]"),
  make_option(c("--width"), type="double", default=10, help="Plot width in inches [default: %default]"),
  make_option(c("--height"), type="double", default=7,  help="Plot height in inches [default: %default]"),
  make_option(c("--max_labels"), type="integer", default=3,
            help="Max contig labels per tile; adds '+' if truncated [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---------- Helpers ----------
trim_ws <- function(x) {
  x %>% str_replace_all("\\s+", " ") %>% str_trim()
}

# Clean gene names a bit: drop trailing allele suffixes like _1 or /variants if present
clean_gene <- function(x) {
  x %>% str_replace("_\\d+$","") %>% trim_ws()
}

read_abricate <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE) %>%
    clean_names()
  # Ensure required columns exist
  required <- c("sample_id","database","product","resistance","gene")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Missing required columns in input: ", paste(missing, collapse=", "))
  }
  df
}

# pick the contig/sequence column name from common variants
pick_seq_col <- function(df) {
  cands <- c("sequence","contig","seqid","contig_id","seq_name","scaffold")
  hit <- intersect(cands, names(df))
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

# extract a friendlier label (numbers if present, else the raw token)
make_contig_label <- function(x, max_n = 3) {
  u <- unique(na.omit(x))
  if (length(u) == 0) return("")
  nums <- stringr::str_extract(u, "\\d+")
  nums <- ifelse(is.na(nums), u, nums)
  nums <- as.character(nums)
  if (length(nums) > max_n) {
    paste0(paste(nums[1:max_n], collapse=","), "+")
  } else {
    paste(nums, collapse=",")
  }
}

# ---------- Load data ----------
dat <- read_abricate(opt$input)

# Tag types
dat <- dat %>%
  mutate(
    record_type = case_when(
      database %in% c("resfinder","card","ncbi","argannot") ~ "AMR",
      database %in% c("plasmidfinder") ~ "PLASMID",
      TRUE ~ toupper(database)
    )
  )

seq_col <- pick_seq_col(dat)
if (is.na(seq_col)) {
  message("No contig/sequence column detected; heatmap tiles will have no contig labels.")
}

# ---------- 1) Distribution: antibiotics (from RESISTANCE) & plasmid types (from PRODUCT) ----------
# Antibiotics: split semicolon list in RESISTANCE for AMR records
abx_long <- dat %>%
  filter(record_type == "AMR", !is.na(resistance), resistance != "") %>%
  separate_rows(resistance, sep = ";") %>%
  mutate(resistance = trim_ws(resistance)) %>%
  filter(resistance != "")

p_abx_dist <- abx_long %>%
  count(resistance, sort = TRUE) %>%
  mutate(resistance = fct_reorder(resistance, n)) %>%
  ggplot(aes(x = resistance, y = n)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Distribution of Antibiotic Resistance Annotations (from RESISTANCE)",
    x = "Antibiotic / Drug class", y = "Count",
    caption = "Antibiotic resistance profiles extracted from RESISTANCE field in ResFinder database"  # Comment out if not needed
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(opt$outdir, "01_distribution_antibiotic_resistance.png"),
       p_abx_dist, width = opt$width, height = opt$height, dpi = 300)

# Plasmid types: PRODUCT field for plasmidfinder
plasmid_long <- dat %>%
  filter(record_type == "PLASMID", !is.na(product), product != "") %>%
  mutate(plasmid_type = trim_ws(product))

p_plasmid_dist <- plasmid_long %>%
  count(plasmid_type, sort = TRUE) %>%
  mutate(plasmid_type = fct_reorder(plasmid_type, n)) %>%
  ggplot(aes(x = plasmid_type, y = n)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Distribution of Plasmid Replicon Types (PlasmidFinder PRODUCT)",
    x = "Replicon / Plasmid type", y = "Count",
    caption = "Plasmid replicon types identified by PlasmidFinder database"  # Comment out if not needed
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(opt$outdir, "02_distribution_plasmid_types.png"),
       p_plasmid_dist, width = opt$width, height = opt$height, dpi = 300)

# ---------- 2) Distribution of AMR genes ----------
gene_field <- if (tolower(opt$gene_column) == "gene") "gene" else "product"

amr_genes <- dat %>%
  filter(record_type == "AMR", !is.na(.data[[gene_field]]), .data[[gene_field]] != "") %>%
  transmute(amr_gene = .data[[gene_field]] %>% clean_gene())

top_genes <- amr_genes %>%
  count(amr_gene, sort = TRUE) %>%
  slice_head(n = opt$topN)

p_amr_dist <- top_genes %>%
  mutate(amr_gene = fct_reorder(amr_gene, n)) %>%
  ggplot(aes(x = amr_gene, y = n)) +
  geom_col() +
  coord_flip() +
  labs(
    title = paste0("Top ", nrow(top_genes), " AMR Genes (", toupper(gene_field), ")"),
    x = "AMR gene", y = "Count",
    caption = paste0("Most prevalent AMR genes detected across all samples (field: ", toupper(gene_field), ")")  # Comment out if not needed
  ) +
  theme_minimal(base_size = 12)

ggsave(file.path(opt$outdir, "03_distribution_amr_genes_topN.png"),
       p_amr_dist, width = opt$width, height = opt$height, dpi = 300)

# ---------- 3) Create unified sample ordering ----------
# Calculate richness for AMR genes
amr_richness <- dat %>%
  filter(record_type == "AMR", !is.na(.data[[gene_field]]), .data[[gene_field]] != "") %>%
  transmute(sample_id, amr_gene = clean_gene(.data[[gene_field]])) %>%
  distinct() %>%
  count(sample_id, name = "amr_richness")

# Calculate richness for plasmids
plasmid_richness <- dat %>%
  filter(record_type == "PLASMID", !is.na(product), product != "") %>%
  transmute(sample_id, plasmid_type = trim_ws(product)) %>%
  distinct() %>%
  count(sample_id, name = "plasmid_richness")

# Combine and create unified sample order
unified_sample_order <- full_join(amr_richness, plasmid_richness, by = "sample_id") %>%
  mutate(
    amr_richness = coalesce(amr_richness, 0L),
    plasmid_richness = coalesce(plasmid_richness, 0L),
    total_richness = amr_richness + plasmid_richness
  ) %>%
  arrange(desc(total_richness), desc(amr_richness), desc(plasmid_richness)) %>%
  pull(sample_id)

# ---------- 4) AMR gene heatmap with contig annotations ----------
seq_col_explicit <- if ("sequence" %in% names(dat)) "sequence" else pick_seq_col(dat)

amr_hits <- dat %>%
  filter(record_type == "AMR",
         !is.na(.data[[gene_field]]), .data[[gene_field]] != "") %>%
  transmute(
    sample_id,
    amr_gene = clean_gene(.data[[gene_field]]),
    seq_lbl  = if (!is.na(seq_col_explicit)) .data[[seq_col_explicit]] else NA_character_
  ) %>%
  distinct()

if (nrow(amr_hits) > 0) {
  # Presence/absence matrix
  amr_pa <- amr_hits %>%
    transmute(sample_id, amr_gene, present = 1L) %>%
    distinct()

  # Ordering (use unified sample order)
  gene_order_amr <- amr_pa %>% count(amr_gene) %>%
    arrange(desc(n)) %>% pull(amr_gene)

  # Build contig labels
  amr_labels <- amr_hits %>%
    group_by(sample_id, amr_gene) %>%
    summarise(contig_label = make_contig_label(seq_lbl, max_n = opt$max_labels),
              .groups = "drop") %>%
    mutate(contig_label = ifelse(contig_label == "", "·", contig_label))

  # Debug output
  try(readr::write_tsv(amr_labels %>% arrange(sample_id, amr_gene),
                       file.path(opt$outdir, "debug_amr_contig_labels.tsv")),
      silent = TRUE)

  # Join and plot
  amr_heat_df <- amr_pa %>%
    left_join(amr_labels, by = c("sample_id","amr_gene")) %>%
    mutate(
      sample_id = factor(sample_id, levels = unified_sample_order),
      amr_gene  = factor(amr_gene,  levels = gene_order_amr)
    )

  p_amr_heat <- ggplot(amr_heat_df, aes(x = sample_id, y = amr_gene, fill = present)) +
    geom_tile(alpha = 0.9) +
    geom_text(aes(label = contig_label),
              color = "white", size = 4, fontface = "bold", na.rm = TRUE) +
    scale_fill_gradient(limits = c(0,1), low = "grey80", high = "steelblue", guide = "none") +
    labs(title = "AMR Gene Presence per Sample (labels = contig IDs/nums)",
         x = "Sample", y = "AMR gene",
         caption = "White labels indicate contig numbers where AMR genes were detected. Samples ordered by total feature richness."  # Comment out if not needed
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  ggsave(file.path(opt$outdir, "04_heatmap_amr_genes_per_sample.png"),
         p_amr_heat, width = 10, height = 10, dpi = 300)
}

# ---------- 4) Plasmid replicon heatmap with contig annotations ----------
plasmid_hits <- dat %>%
  filter(record_type == "PLASMID",
         !is.na(product), product != "") %>%
  transmute(
    sample_id,
    plasmid_type = trim_ws(product),
    seq_lbl      = if (!is.na(seq_col)) .data[[seq_col]] else NA_character_
  ) %>%
  distinct()

if (nrow(plasmid_hits) > 0) {
  plasmid_pa <- plasmid_hits %>%
    transmute(sample_id, plasmid_type, present = 1L) %>%
    distinct()

  # Use unified sample order for plasmids too
  type_order_pl   <- plasmid_pa %>% count(plasmid_type) %>% 
    arrange(desc(n)) %>% pull(plasmid_type)

  plasmid_labels <- plasmid_hits %>%
    group_by(sample_id, plasmid_type) %>%
    summarise(contig_label = make_contig_label(seq_lbl, max_n = opt$max_labels),
              .groups = "drop") %>%
    mutate(contig_label = ifelse(contig_label == "", "·", contig_label))

  plasmid_heat_df <- plasmid_pa %>%
    left_join(plasmid_labels, by = c("sample_id","plasmid_type")) %>%
    mutate(
      sample_id    = factor(sample_id, levels = unified_sample_order),
      plasmid_type = factor(plasmid_type, levels = type_order_pl)
    )

  p_plasmid_heat <- ggplot(plasmid_heat_df, aes(x = sample_id, y = plasmid_type, fill = present)) +
    geom_tile(alpha = 0.9) +
    geom_text(aes(label = contig_label),
              color = "white", size = 4, fontface = "bold", na.rm = TRUE) +
    scale_fill_gradient(limits = c(0,1), low = "grey80", high = "steelblue", guide = "none") +
    labs(title = "Plasmid Replicon Types per Sample (labels = contig IDs/nums)",
         x = "Sample", y = "Replicon / type",
         caption = "White labels indicate contig numbers where plasmid replicons were detected. Sample order matches AMR heatmap."  # Comment out if not needed
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  ggsave(file.path(opt$outdir, "05_heatmap_plasmid_types_per_sample.png"),
         p_plasmid_heat, width = 10, height = 10, dpi = 300)
}

# ---------- Save simple TSV summaries ----------
# Antibiotic counts
if (nrow(abx_long) > 0) {
  abx_long %>% count(resistance, sort = TRUE) %>%
    write_tsv(file.path(opt$outdir, "summary_antibiotic_resistance_counts.tsv"))
}

# Plasmid type counts
if (nrow(plasmid_long) > 0) {
  plasmid_long %>% count(plasmid_type, sort = TRUE) %>%
    write_tsv(file.path(opt$outdir, "summary_plasmid_type_counts.tsv"))
}

# AMR gene counts
amr_genes %>% count(amr_gene, sort = TRUE) %>%
  write_tsv(file.path(opt$outdir, "summary_amr_gene_counts.tsv"))

# ---------- 5) Per-sample summary table (lists + counts) ----------
amr_by_sample <- dat %>%
  filter(record_type == "AMR",
         !is.na(.data[[gene_field]]), .data[[gene_field]] != "") %>%
  transmute(sample_id, amr_gene = clean_gene(.data[[gene_field]])) %>%
  distinct() %>%
  group_by(sample_id) %>%
  summarise(
    amr_genes        = paste(sort(amr_gene), collapse = ";"),
    amr_gene_count   = n(),
    .groups = "drop"
  )

plasmid_by_sample <- dat %>%
  filter(record_type == "PLASMID",
         !is.na(product), product != "") %>%
  transmute(sample_id, plasmid_replicon = trim_ws(product)) %>%
  distinct() %>%
  group_by(sample_id) %>%
  summarise(
    plasmid_replicons      = paste(sort(plasmid_replicon), collapse = ";"),
    plasmid_replicon_count = n(),
    .groups = "drop"
  )

per_sample_summary <- full_join(amr_by_sample, plasmid_by_sample, by = "sample_id") %>%
  mutate(
    amr_genes               = coalesce(amr_genes, ""),
    plasmid_replicons       = coalesce(plasmid_replicons, ""),
    amr_gene_count          = coalesce(amr_gene_count, 0L),
    plasmid_replicon_count  = coalesce(plasmid_replicon_count, 0L),
    total_feature_count     = amr_gene_count + plasmid_replicon_count
  ) %>%
  arrange(desc(total_feature_count), sample_id)

write_tsv(per_sample_summary,
          file.path(opt$outdir, "per_sample_amr_plasmid_summary.tsv"))

message("Done. Plots and summaries written to: ", opt$outdir)