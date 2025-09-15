#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) {
    install.packages("optparse", repos = "https://cloud.r-project.org")
  }
  if (!requireNamespace("tidyverse", quietly = TRUE)) {
    install.packages("tidyverse", repos = "https://cloud.r-project.org")
  }
  library(optparse)
  library(tidyverse)  # dplyr, readr, ggplot2, forcats, stringr, tibble, purrr
})

# ---------------------------
# CLI
# ---------------------------
option_list <- list(
  make_option(c("-i","--input"),  type="character", default="biosample.csv",
              help="Input CSV file [default: %default]"),
  make_option(c("-o","--outdir"), type="character", default="plots",
              help="Output directory [default: %default]"),
  make_option(c("-s","--source"), type="character", default="urine",
              help="Comma-separated source terms to subset (e.g., 'urine,blood') [default: %default]")
)
opt <- parse_args(OptionParser(option_list=option_list))

in_csv  <- opt$input
out_dir <- opt$outdir
source_filter_raw <- opt$source

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Ensure parsed terms exist and are character vector
source_terms <- stringr::str_split(source_filter_raw, ",", simplify = TRUE) |>
  as.character() |> trimws()
source_terms <- source_terms[source_terms != ""]
if (length(source_terms) == 0) {
  source_terms <- c("urine")
}

# Containers
subsets_by_term <- list()
term_hits_summaries <- list()

# Helper: SI labels compatible with modern 'scales'
if (!requireNamespace("scales", quietly = TRUE)) {
  install.packages("scales", repos = "https://cloud.r-project.org")
}
si_labels <- scales::label_number(scale_cut = scales::cut_short_scale())

# ---------------------------
# Helpers
# ---------------------------

# Coalesce a column from a set of name candidates (regex OK, case-insensitive)
pick_col <- function(df, candidates) {
  nms <- names(df)
  for (pat in candidates) {
    idx <- which(stringr::str_detect(nms, regex(pat, ignore_case = TRUE)))
    if (length(idx) > 0) return(nms[idx[1]])
  }
  return(NA_character_)
}

# Make a clean DF with standardized columns if present
standardize_biosample <- function(df) {
  nm <- names(df)
  get_col <- function(patterns) {
    col <- pick_col(df, patterns)
    if (is.na(col)) return(rep(NA, nrow(df)))
    df[[col]]
  }
  tibble::tibble(
    accession        = get_col(c("^accession$","^acc$","^biosample[_ ]?id$","^sample[_ ]?accession$")),
    title            = get_col(c("^title$","^sample[_ ]?title$")),
    organism         = get_col(c("^organism$","^taxon","^species$")),
    host             = get_col(c("^host$","host.*name","host.*species")),
    sample_source    = get_col(c("^sample[_ ]?source$","^isolation[_ ]?source$","^source$","specimen.*type","^material$","body[_ ]?site","anatomical")),
    country          = get_col(c("^country$","geographic.*country","location","^nation$")),
    collection_year  = suppressWarnings(as.integer(get_col(c("^collection[_ ]?year$","collection.*date","year"))))
  )
}

# Read input
message("Reading: ", in_csv)
raw <- suppressMessages(readr::read_csv(in_csv, show_col_types = FALSE))

# Expand df_exp with standardized columns, but keep all original columns as well
df_exp <- raw
std <- standardize_biosample(raw)
for (col in names(std)) {
  if (!col %in% names(df_exp)) df_exp[[col]] <- std[[col]]
}

# Isolation explicit column (if a more exact source column exists)
col_iso <- pick_col(raw, c("^isolation[_ ]?source$","isolation source","sample_source_detail",
                           "specimen_source","sample_type","specimen","material",
                           "body[_ ]?site","anatomical"))
df_exp <- df_exp |>
  mutate(
    isolation_source_explicit = if (!is.na(col_iso)) raw[[col_iso]] else sample_source,
    isolation_source_explicit = stringr::str_squish(as.character(isolation_source_explicit))
  )

# Searchable text columns
search_cols <- c("sample_source","isolation_source_explicit","isolation_source","sample_source_detail","specimen_source","sample_type","specimen","material","body_site","anatomical")
search_cols <- intersect(search_cols, names(df_exp))

# Also auto-include any column names containing "source" or "specimen"
auto_cols <- names(df_exp)[stringr::str_detect(names(df_exp), regex("(source|specimen)", ignore_case = TRUE))]
search_cols <- unique(c(search_cols, auto_cols))

if (length(search_cols) == 0) {
  message("No obvious search columns found; considering all character columns.")
  search_cols <- names(df_exp)[vapply(df_exp, inherits, logical(1), what = "character")]
}

# ---------------------------
# Regex plan
# ---------------------------
canon_terms <- unique(tolower(source_terms))

# For blood: user wants the broadest match: any substring containing 'blood' (case-insensitive)
# For others: use the term itself as substring (e.g., 'urine'), case-insensitive
per_term_regex <- setNames(vector("character", length(canon_terms)), canon_terms)
for (i in seq_along(canon_terms)) {
  tm <- canon_terms[[i]]
  if (tm == "blood") {
    per_term_regex[[tm]] <- "(?i)blood"   # simplest, broadest
  } else if (tm == "urine") {
    per_term_regex[[tm]] <- "(?i)urine|(?i)urinary|(?i)uti"  # a couple of obvious variants
  } else {
    # generic fallback: the literal term as substring, case-insensitive
    esc <- stringr::str_replace_all(tm, "([\\^$.|()?*+\\[\\]{}\\\\])", "\\\\\\1")
    per_term_regex[[tm]] <- paste0("(?i)", esc)
  }
}

# Combined regex across all terms
combined_regex <- paste(per_term_regex, collapse = "|")

# ---------------------------
# Build subsets
# ---------------------------

# Combined subset (any term in any search column)
combined_subset <- df_exp %>%
  mutate(across(all_of(search_cols), ~ tolower(as.character(.x)))) %>%
  filter(if_any(all_of(search_cols),
                ~ stringr::str_detect(coalesce(.x, ""), regex(combined_regex)) )) %>%
  distinct(accession, .keep_all = TRUE)

subset_tag <- tolower(paste(canon_terms, collapse = "_")) |> gsub("\\W+","_", x = _)
out_csv_combined <- file.path(out_dir, paste0("samples_", subset_tag, "_only.csv"))
if (nrow(combined_subset) > 0) {
  readr::write_csv(combined_subset, out_csv_combined)
  message("Saved combined subset (", paste(canon_terms, collapse=","), "): ", out_csv_combined)
} else {
  message("No rows matched any of: ", paste(canon_terms, collapse=", "))
}

# Per-term subsets
for (term in canon_terms) {
  term_regex <- per_term_regex[[term]]
  subset_term <- df_exp %>%
    mutate(across(all_of(search_cols), ~ tolower(as.character(.x)))) %>%
    filter(if_any(all_of(search_cols),
                  ~ stringr::str_detect(coalesce(.x, ""), regex(term_regex)) )) %>%
    distinct(accession, .keep_all = TRUE)

  term_tag <- tolower(gsub("\\W+","_", term))
  subsets_by_term[[term_tag]] <- subset_term

  # Diagnostics for this term
  hits_df <- tibble(
    term = term,
    column = search_cols,
    n_matches = vapply(search_cols, function(cc) {
      vv <- tolower(as.character(df_exp[[cc]]))
      sum(stringr::str_detect(coalesce(vv, ""), regex(term_regex)))
    }, integer(1))
  )
  term_hits_summaries[[term_tag]] <- hits_df

  if (nrow(subset_term) > 0) {
    out_csv_t <- file.path(out_dir, paste0("samples_", term_tag, "_only.csv"))
    readr::write_csv(subset_term, out_csv_t)
    message("Saved ", term, "-only subset: ", out_csv_t)

    # ---- Plots for this term ----
    # Organism counts
    if ("organism" %in% names(subset_term)) {
      tab_org <- subset_term %>%
        filter(!is.na(organism), organism != "") %>%
        count(organism, sort = TRUE) %>%
        slice_head(n = 20)
      if (nrow(tab_org) > 0) {
        p1 <- tab_org %>%
          mutate(organism = forcats::fct_reorder(organism, n)) %>%
          ggplot(aes(organism, n)) +
          geom_col() + coord_flip() +
          labs(title = paste("Subset:", term, "— Organisms (top 20)"),
               x = "Organism", y = "Count") +
          scale_y_continuous(labels = si_labels) +
          theme_minimal(base_size = 12)
        ggsave(file.path(out_dir, paste0("SUB1_counts_organism_", term_tag, "_only.png")),
               p1, width = 8, height = 6, dpi = 300)
      }
    }

    # Host counts
    if ("host" %in% names(subset_term)) {
      tab_host <- subset_term %>%
        filter(!is.na(host), host != "") %>%
        count(host, sort = TRUE) %>%
        slice_head(n = 20)
      if (nrow(tab_host) > 0) {
        p2 <- tab_host %>%
          mutate(host = forcats::fct_reorder(host, n)) %>%
          ggplot(aes(host, n)) +
          geom_col() + coord_flip() +
          labs(title = paste("Subset:", term, "— Hosts (top 20)"),
               x = "Host", y = "Count") +
          scale_y_continuous(labels = si_labels) +
          theme_minimal(base_size = 12)
        ggsave(file.path(out_dir, paste0("SUB2_counts_host_", term_tag, "_only.png")),
               p2, width = 8, height = 6, dpi = 300)
      }
    }

    # Sample source labels
    if ("sample_source" %in% names(subset_term)) {
      tab_src <- subset_term %>%
        filter(!is.na(sample_source), sample_source != "") %>%
        count(sample_source, sort = TRUE) %>%
        slice_head(n = 20)
      if (nrow(tab_src) > 0) {
        p3 <- tab_src %>%
          mutate(sample_source = forcats::fct_reorder(sample_source, n)) %>%
          ggplot(aes(sample_source, n)) +
          geom_col() + coord_flip() +
          labs(title = paste("Subset:", term, "— Sample source labels (top 20)"),
               x = "Sample source", y = "Count") +
          scale_y_continuous(labels = si_labels) +
          theme_minimal(base_size = 12)
        ggsave(file.path(out_dir, paste0("SUB3_counts_sample_source_", term_tag, "_only.png")),
               p3, width = 8, height = 6, dpi = 300)
      }
    }

    # Collection year
    if ("collection_year" %in% names(subset_term)) {
      p4 <- subset_term %>%
        filter(!is.na(collection_year), collection_year > 1900, collection_year < 2100) %>%
        ggplot(aes(x = collection_year)) +
        geom_histogram(binwidth = 1, boundary = 0, closed = "right") +
        scale_x_continuous(breaks = scales::pretty_breaks()) +
        labs(title = paste("Subset:", term, "— Collection year"),
             x = "Year", y = "Samples") +
        theme_minimal(base_size = 12)
      ggsave(file.path(out_dir, paste0("SUB4_hist_year_", term_tag, "_only.png")),
             p4, width = 8, height = 6, dpi = 300)
    }

    # Organism by Host
    if (all(c("organism","host") %in% names(subset_term))) {
      tab <- subset_term %>%
        filter(!is.na(organism), organism != "", !is.na(host), host != "") %>%
        count(organism, host, sort = TRUE)
      if (nrow(tab) > 0) {
        top_orgs <- tab %>% group_by(organism) %>% summarise(n = sum(n), .groups="drop") %>%
          slice_max(n, n = 15) %>% pull(organism)
        p5 <- tab %>%
          mutate(organism = if_else(organism %in% top_orgs, organism, "Other")) %>%
          group_by(organism, host) %>% summarise(n = sum(n), .groups="drop") %>%
          mutate(organism = forcats::fct_reorder(organism, n, sum)) %>%
          ggplot(aes(x = organism, y = n, fill = host)) +
          geom_col(position = "stack") + coord_flip() +
          labs(title = paste("Subset:", term, "— Organism × Host"),
               x = "Organism", y = "Count", fill = "Host") +
          theme_minimal(base_size = 12)
        ggsave(file.path(out_dir, paste0("SUB5_stack_organism_by_host_", term_tag, "_only.png")),
               p5, width = 9, height = 7, dpi = 300)
      }
    }

    # Location / country counts
    if ("country" %in% names(subset_term)) {
      tab_loc <- subset_term %>%
        filter(!is.na(country), country != "") %>%
        count(country, sort = TRUE) %>%
        slice_head(n = 30)
      if (nrow(tab_loc) > 0) {
        p6 <- tab_loc %>%
          mutate(country = as.character(country)) %>%
          mutate(country = forcats::fct_reorder(country, n)) %>%
          ggplot(aes(country, n)) +
          geom_col() + coord_flip() +
          labs(title = paste("Subset:", term, "— Samples by location"),
               x = "Location / Country", y = "Count") +
          theme_minimal(base_size = 12)
        ggsave(file.path(out_dir, paste0("SUB6_counts_location_", term_tag, "_only.png")),
               p6, width = 8, height = 6, dpi = 300)
      }
    }

  } else {
    message("No samples matched source term: ", term)
  }
}


# ---------------------------
# Combined comparison plots: any N terms (>= 2)
# ---------------------------
avail_terms <- names(subsets_by_term)
have_any <- avail_terms[vapply(subsets_by_term[avail_terms], function(df) !is.null(df) && nrow(df) > 0, logical(1))]

if (length(have_any) >= 2) {
  combined_df_all <- dplyr::bind_rows(lapply(have_any, function(tg) {
    df <- subsets_by_term[[tg]]
    if (!is.null(df) && nrow(df) > 0) { df$source_term <- tg; df } else { NULL }
  }))

  if (!is.null(combined_df_all) && nrow(combined_df_all) > 0) {

    # Ensure country is character if present
    if ("country" %in% names(combined_df_all)) {
      combined_df_all <- dplyr::mutate(combined_df_all, country = as.character(country))
    }

    # 1) Organisms by source_term (top 15 overall)
    if (all(c("organism","source_term") %in% names(combined_df_all))) {
      tab_org_all <- combined_df_all %>%
        dplyr::filter(!is.na(organism), organism != "") %>%
        dplyr::count(source_term, organism, name = "n", sort = TRUE)

      top_orgs_all <- tab_org_all %>%
        dplyr::group_by(organism) %>% dplyr::summarise(N = sum(n), .groups = "drop") %>%
        dplyr::slice_max(N, n = 15) %>% dplyr::pull(organism)

      p_comb_org_all <- tab_org_all %>%
        dplyr::mutate(organism = dplyr::if_else(organism %in% top_orgs_all, organism, "Other")) %>%
        dplyr::group_by(source_term, organism) %>% dplyr::summarise(n = sum(n), .groups = "drop") %>%
        dplyr::mutate(organism = forcats::fct_reorder(organism, n, sum)) %>%
        ggplot2::ggplot(ggplot2::aes(x = organism, y = n, fill = source_term)) +
        ggplot2::geom_col(position = "dodge") +
        ggplot2::coord_flip() +
        ggplot2::labs(title = "Organisms by source term",
                      x = "Organism", y = "Count", fill = "Source term") +
        ggplot2::scale_y_continuous(labels = si_labels) +
        ggplot2::theme_minimal(base_size = 12)
      ggplot2::ggsave(file.path(out_dir, "COMBO_counts_organism_all_terms.png"),
                      p_comb_org_all, width = 10, height = 8, dpi = 300)
    }

    # 2) Collection year histograms faceted by source_term
    if (all(c("collection_year","source_term") %in% names(combined_df_all))) {
      p_comb_year_all <- combined_df_all %>%
        dplyr::filter(!is.na(collection_year), collection_year > 1900, collection_year < 2100) %>%
        ggplot2::ggplot(ggplot2::aes(x = collection_year)) +
        ggplot2::geom_histogram(binwidth = 1, boundary = 0, closed = "right") +
        ggplot2::facet_wrap(~ source_term, ncol = 1, scales = "free_y") +
        ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
        ggplot2::labs(title = "Collection year by source term",
                      x = "Year", y = "Samples") +
        ggplot2::theme_minimal(base_size = 12)
      ggplot2::ggsave(file.path(out_dir, "COMBO_hist_year_all_terms.png"),
                      p_comb_year_all, width = 10, height = 10, dpi = 300)
    }

    # 3) Location / country counts by source_term (top 20 overall)
    if (all(c("country","source_term") %in% names(combined_df_all))) {
      tab_loc_all <- combined_df_all %>%
        dplyr::mutate(country = as.character(country)) %>%
        dplyr::filter(!is.na(country), country != "") %>%
        dplyr::count(source_term, country, name = "n", sort = TRUE)

      top_cty_all <- tab_loc_all %>%
        dplyr::group_by(country) %>% dplyr::summarise(N = sum(n), .groups = "drop") %>%
        dplyr::slice_max(N, n = 20) %>% dplyr::pull(country)

      p_comb_loc_all <- tab_loc_all %>%
        dplyr::mutate(country = as.character(country)) %>%
        dplyr::mutate(country = dplyr::if_else(country %in% top_cty_all, country, "Other")) %>%
        dplyr::group_by(source_term, country) %>% dplyr::summarise(n = sum(n), .groups = "drop") %>%
        dplyr::mutate(country = forcats::fct_reorder(country, n, sum)) %>%
        ggplot2::ggplot(ggplot2::aes(x = country, y = n, fill = source_term)) +
        ggplot2::geom_col(position = "dodge") +
        ggplot2::coord_flip() +
        ggplot2::labs(title = "Samples by location (top 20), by source term",
                      x = "Location / Country", y = "Count", fill = "Source term") +
        ggplot2::theme_minimal(base_size = 12)
      ggplot2::ggsave(file.path(out_dir, "COMBO_counts_location_all_terms.png"),
                      p_comb_loc_all, width = 11, height = 9, dpi = 300)
    }
  }
} else {
  message("Generic combined comparison skipped (need >= 2 non-empty terms).")
}

# Backward-compat: if exactly two terms are urine+blood, also write the legacy filenames
if (setequal(have_any, c("urine","blood"))) {
  if (file.exists(file.path(out_dir, "COMBO_counts_organism_all_terms.png"))) {
    file.copy(file.path(out_dir, "COMBO_counts_organism_all_terms.png"),
              file.path(out_dir, "COMBO_counts_organism_urine_vs_blood.png"), overwrite = TRUE)
  }
  if (file.exists(file.path(out_dir, "COMBO_hist_year_all_terms.png"))) {
    file.copy(file.path(out_dir, "COMBO_hist_year_all_terms.png"),
              file.path(out_dir, "COMBO_hist_year_urine_vs_blood.png"), overwrite = TRUE)
  }
  if (file.exists(file.path(out_dir, "COMBO_counts_location_all_terms.png"))) {
    file.copy(file.path(out_dir, "COMBO_counts_location_all_terms.png"),
              file.path(out_dir, "COMBO_counts_location_urine_vs_blood.png"), overwrite = TRUE)
  }
}
# ---------------------------
# Diagnostics summary
# ---------------------------
if (length(term_hits_summaries) > 0) {
  diag_df <- dplyr::bind_rows(term_hits_summaries, .id = "term_tag")
  if (!is.null(diag_df) && nrow(diag_df) > 0) {
    readr::write_csv(diag_df, file.path(out_dir, "term_hits_summary.csv"))
    message("Wrote diagnostics: ", file.path(out_dir, "term_hits_summary.csv"))
    print(diag_df %>% group_by(term) %>% summarise(total_matches = sum(n_matches), .groups="drop"))
  }
}
