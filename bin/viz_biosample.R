# ---- Setup ----
# Install packages once if needed:
# install.packages(c("tidyverse","janitor","lubridate","scales"))

library(tidyverse)
library(janitor)
library(lubridate)
library(scales)

# Vector-safe coalesce for character vectors
vcoalesce <- function(x, fill = "") {
  if (is.null(x) || length(x) == 0) return(fill)
  x <- as.character(x)
  x[is.na(x)] <- fill
  x
}

# ====== USER INPUT ======
in_csv <- "biosample.csv"   # <-- change if your file name differs
out_dir <- "plots"          # all figures saved here
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Helpers ----
# Try to pick a column from multiple candidate names (case-insensitive, flexible)
pick_col <- function(df, candidates) {
  # candidates: character vector of regex patterns (case-insensitive)
  cn <- names(df)
  for (pat in candidates) {
    hit <- cn[str_detect(tolower(cn), tolower(pat))]
    if (length(hit) > 0) return(hit[1])
  }
  return(NA_character_)
}

# Split multi-values separated by " | " and unnest
split_multi <- function(x) {
  x %>%
    mutate(across(everything(), ~na_if(.x, ""))) %>%
    separate_rows(value, sep = "\\s*\\|\\s*", convert = FALSE)
}

# ---- Read & normalize ----
raw <- readr::read_csv(in_csv, show_col_types = FALSE) %>% clean_names()

# Identify likely columns
col_accession <- pick_col(raw, c("^accession$","^sam(n|e|d)\\w*","^biosample"))
col_title     <- pick_col(raw, c("^title$"))
col_org       <- pick_col(raw, c("^organism$","^species$"))
col_host      <- pick_col(raw, c("^host$"))
col_source    <- pick_col(raw, c("^samplesource$","^sample_source$",
                                 "^isolation_source$","(^|_)source$","^sample_type$","^isolation$"))
col_loc       <- pick_col(raw, c("^geographic.*", "^country$", "^location$", "geo", "isolation_country"))
col_year      <- pick_col(raw, c("^collectionyear$","^collection_year$","^year$","^collection_date$","^collectiondate$"))

message("Detected columns:\n",
        "  Accession:      ", col_accession %||% "—", "\n",
        "  Title:          ", col_title %||% "—", "\n",
        "  Organism:       ", col_org %||% "—", "\n",
        "  Host:           ", col_host %||% "—", "\n",
        "  Sample Source:  ", col_source %||% "—", "\n",
        "  Location:       ", col_loc %||% "—", "\n",
        "  CollectionYear: ", col_year %||% "—", "\n")

`%||%` <- function(a,b) if (is.null(a) || is.na(a)) b else a

df <- raw %>%
  mutate(
    accession      = if (!is.na(col_accession)) .data[[col_accession]] else NA_character_,
    title          = if (!is.na(col_title)) .data[[col_title]] else NA_character_,
    organism       = if (!is.na(col_org)) .data[[col_org]] else NA_character_,
    host           = if (!is.na(col_host)) .data[[col_host]] else NA_character_,
    sample_source  = if (!is.na(col_source)) .data[[col_source]] else NA_character_,
    location_raw   = if (!is.na(col_loc)) .data[[col_loc]] else NA_character_,
    collection_raw = if (!is.na(col_year)) .data[[col_year]] else NA_character_
  )

# ---- Derive/clean fields ----
# 1) Collection year: if a full date present, extract year; else keep numeric year if present.
extract_year <- function(x) {
  if (is.null(x)) return(NA_integer_)
  x_chr <- as.character(x)
  # Try parse date first
  y <- suppressWarnings(year(ymd(x_chr)))
  # If parsing failed, try regex for a 4-digit year (19xx/20xx)
  y2 <- ifelse(is.na(y),
               str_extract(x_chr, "(19|20)\\d{2}"),
               as.character(y))
  as.integer(y2)
}

df <- df %>%
  mutate(
    collection_year = extract_year(collection_raw),
    organism = str_squish(as.character(organism)),
    host = str_squish(as.character(host)),
    sample_source = str_squish(as.character(sample_source)),
    location_raw = str_squish(as.character(location_raw))
  )

# For multi-valued fields separated by " | ", we’ll also create exploded versions
explode_multi <- function(df, col) {
  if (is.na(col)) return(df %>% mutate("{col}" := NA_character_))
  d <- df %>%
    mutate(value = .data[[col]]) %>%
    mutate(value = ifelse(is.na(value) | value == "", NA_character_, value)) %>%
    separate_rows(value, sep = "\\s*\\|\\s*") %>%
    mutate("{col}" := value) %>%
    select(-value)
  return(d)
}

df_org  <- explode_multi(df, "organism")
df_host <- explode_multi(df_org, "host")
df_src  <- explode_multi(df_host, "sample_source")
df_loc  <- explode_multi(df_src, "location_raw")
df_exp  <- df_loc

# Make nicer labels
nice_counts_plot <- function(data, aes_x, top_n = 20, title = "", x_lab = "",
                             y_lab = "Count", filename) {
  tab <- data %>%
    filter(!is.na({{aes_x}}), {{aes_x}} != "") %>%
    count({{aes_x}}, sort = TRUE) %>%
    slice_head(n = top_n)
  if (nrow(tab) == 0) {
    message("Skipping plot (no data): ", title)
    return(invisible(NULL))
  }
  p <- tab %>%
    mutate({{aes_x}} := fct_reorder({{aes_x}}, n)) %>%
    ggplot(aes({{aes_x}}, n)) +
    geom_col() +
    coord_flip() +
    labs(title = title, x = x_lab, y = y_lab) +
    scale_y_continuous(
      labels = label_number(scale_cut = cut_si(" "))
    ) +
    theme_minimal(base_size = 12)
  ggsave(file.path(out_dir, filename), p, width = 8, height = 6, dpi = 300)
  print(p)
}

# ---- 1) Organism counts (top 20) ----
if (!is.na(col_org)) {
  nice_counts_plot(df_exp, organism, 20,
                   title = "Top organisms",
                   x_lab = "Organism",
                   filename = "01_counts_organism_top20.png")
}

# ---- 2) Sample Source counts ----
if (!is.na(col_source)) {
  nice_counts_plot(df_exp, sample_source, 30,
                   title = "Sample sources",
                   x_lab = "Sample source",
                   filename = "02_counts_sample_source.png")
}

# ---- 3) Host counts ----
if (!is.na(col_host)) {
  nice_counts_plot(df_exp, host, 30,
                   title = "Hosts",
                   x_lab = "Host",
                   filename = "03_counts_host.png")
}

# ---- 4) Collection year histogram ----
if (!is.na(col_year) || "collection_year" %in% names(df_exp)) {
  p_year <- df_exp %>%
    filter(!is.na(collection_year), collection_year > 1900, collection_year < 2100) %>%
    ggplot(aes(x = collection_year)) +
    geom_histogram(binwidth = 1, boundary = 0, closed = "right") +
    scale_x_continuous(breaks = pretty_breaks()) +
    labs(title = "Collection year distribution", x = "Year", y = "Samples") +
    theme_minimal(base_size = 12)
  if (nrow(ggplot_build(p_year)$data[[1]]) > 0) {
    ggsave(file.path(out_dir, "04_hist_collection_year.png"), p_year, width = 8, height = 5, dpi = 300)
    print(p_year)
  } else {
    message("Skipping year histogram (no valid year data).")
  }
}

# ---- 5) Stacked: Organism by Sample Source ----
if (!is.na(col_org) && !is.na(col_source)) {
  tab <- df_exp %>%
    filter(!is.na(organism), organism != "", !is.na(sample_source), sample_source != "") %>%
    count(organism, sample_source, sort = TRUE)
  if (nrow(tab) > 0) {
    top_orgs <- tab %>%
      group_by(organism) %>% summarise(n = sum(n), .groups = "drop") %>%
      slice_max(n, n = 15) %>% pull(organism)
    p_os <- tab %>%
      mutate(organism = if_else(organism %in% top_orgs, organism, "Other")) %>%
      group_by(organism, sample_source) %>% summarise(n = sum(n), .groups = "drop") %>%
      mutate(organism = fct_reorder(organism, n, sum)) %>%
      ggplot(aes(x = organism, y = n, fill = sample_source)) +
      geom_col(position = "stack") +
      coord_flip() +
      labs(title = "Organism × Sample source", x = "Organism", y = "Count", fill = "Sample source") +
      theme_minimal(base_size = 12)
    ggsave(file.path(out_dir, "05_stack_organism_by_source.png"), p_os, width = 9, height = 7, dpi = 300)
    print(p_os)
  }
}

# ---- 6) Stacked: Organism by Host ----
if (!is.na(col_org) && !is.na(col_host)) {
  tab <- df_exp %>%
    filter(!is.na(organism), organism != "", !is.na(host), host != "") %>%
    count(organism, host, sort = TRUE)
  if (nrow(tab) > 0) {
    top_orgs <- tab %>%
      group_by(organism) %>% summarise(n = sum(n), .groups = "drop") %>%
      slice_max(n, n = 15) %>% pull(organism)
    p_oh <- tab %>%
      mutate(organism = if_else(organism %in% top_orgs, organism, "Other")) %>%
      group_by(organism, host) %>% summarise(n = sum(n), .groups = "drop") %>%
      mutate(organism = fct_reorder(organism, n, sum)) %>%
      ggplot(aes(x = organism, y = n, fill = host)) +
      geom_col(position = "stack") +
      coord_flip() +
      labs(title = "Organism × Host", x = "Organism", y = "Count", fill = "Host") +
      theme_minimal(base_size = 12)
    ggsave(file.path(out_dir, "06_stack_organism_by_host.png"), p_oh, width = 9, height = 7, dpi = 300)
    print(p_oh)
  }
}

# ---- 7) Heatmap: Organism × Host ----
if (!is.na(col_org) && !is.na(col_host)) {
  hm <- df_exp %>%
    filter(!is.na(organism), organism != "", !is.na(host), host != "") %>%
    count(organism, host) %>%
    group_by(organism) %>% mutate(total_org = sum(n)) %>% ungroup() %>%
    mutate(organism_top = fct_reorder(organism, total_org))
  if (nrow(hm) > 0) {
    p_hm <- hm %>%
      ggplot(aes(host, organism_top, fill = n)) +
      geom_tile() +
      scale_fill_continuous(type = "viridis") +
      labs(title = "Heatmap: Organism × Host", x = "Host", y = "Organism (by total)", fill = "Count") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(out_dir, "07_heatmap_organism_host.png"), p_hm, width = 10, height = 7, dpi = 300)
    print(p_hm)
  }
}

# ---- 8) Location / Country bar chart  ----
if (!is.na(col_loc)) {
  # Extract a coarse "country-ish" token: take last token after ":" or split by ";" or ","
  country_guess <- df_exp %>%
    mutate(country = location_raw) %>%
    mutate(country = str_replace_all(country, "\\s+", " ")) %>%
    mutate(country = str_replace(country, ".*:\\s*", "")) %>%  # after 'country: X'
    mutate(country = str_split(country, ";|,") %>% map_chr(~ str_squish(.x[[length(.x)]] %||% NA_character_))) %>%
    mutate(country = if_else(country == "", NA_character_, country))

  p_loc <- country_guess %>%
    filter(!is.na(country), country != "") %>%
    count(country, sort = TRUE) %>%
    slice_head(n = 30) %>%
    mutate(country = fct_reorder(country, n)) %>%
    ggplot(aes(country, n)) +
    geom_col() +
    coord_flip() +
    labs(title = "Samples by location", x = "Location / Country ", y = "Count") +
    theme_minimal(base_size = 12)

  if (nrow(ggplot_build(p_loc)$data[[1]]) > 0) {
    ggsave(file.path(out_dir, "08_counts_location.png"), p_loc, width = 8, height = 6, dpi = 300)
    print(p_loc)
  } else {
    message("Skipping location bar chart (no parsable locations).")
  }
}

message("Done. PNGs saved to: ", normalizePath(out_dir))

# ---- Filter: only "urine" samples ----
# pick the most relevant column
col_iso <- pick_col(raw, c("^isolation_source$", "isolation source",
                           "sample_source_detail", "specimen_source",
                           "body[_ ]?site", "anatomical", "isolation"))

df_exp <- df_exp %>%
  mutate(
    isolation_source_explicit = if (!is.na(col_iso)) raw[[col_iso]] else sample_source,
    isolation_source_explicit = str_squish(as.character(isolation_source_explicit)),
    sample_source = str_squish(as.character(sample_source))
  ) %>%
  separate_rows(isolation_source_explicit, sep = "\\s*\\|\\s*") %>%
  separate_rows(sample_source, sep = "\\s*\\|\\s*")

# lowercase versions
df_exp <- df_exp %>%
  mutate(
    iso_lc = tolower(isolation_source_explicit),
    src_lc = tolower(sample_source)
  )

# filter: must mention "urine" in either sample_source or isolation_source
urine_samples <- df_exp %>%
  filter(str_detect(coalesce(src_lc, ""), "\\burine\\b") |
           str_detect(coalesce(iso_lc, ""), "\\burine\\b")) %>%
  distinct(accession, .keep_all = TRUE)

# show + save
print(urine_samples)
if (nrow(urine_samples) > 0) {
  write_csv(urine_samples, file.path(out_dir, "samples_urine_only.csv"))
  message("Saved urine-only subset: ", file.path(out_dir, "samples_urine_only.csv"))
} else {
  message("No urine samples found.")
}

# ---- Filter: only "urine" samples ----
# pick the most relevant column
col_iso <- pick_col(raw, c("^isolation_source$", "isolation source",
                           "sample_source_detail", "specimen_source",
                           "body[_ ]?site", "anatomical", "isolation"))

df_exp <- df_exp %>%
  mutate(
    isolation_source_explicit = if (!is.na(col_iso)) raw[[col_iso]] else sample_source,
    isolation_source_explicit = str_squish(as.character(isolation_source_explicit)),
    sample_source = str_squish(as.character(sample_source))
  ) %>%
  separate_rows(isolation_source_explicit, sep = "\\s*\\|\\s*") %>%
  separate_rows(sample_source, sep = "\\s*\\|\\s*")

# lowercase versions
df_exp <- df_exp %>%
  mutate(
    iso_lc = tolower(isolation_source_explicit),
    src_lc = tolower(sample_source)
  )

# filter: must mention "urine" in either sample_source or isolation_source
urine_samples <- df_exp %>%
  filter(str_detect(coalesce(src_lc, ""), "\\burine\\b") |
           str_detect(coalesce(iso_lc, ""), "\\burine\\b")) %>%
  distinct(accession, .keep_all = TRUE)

# show + save
print(urine_samples)
if (nrow(urine_samples) > 0) {
  write_csv(urine_samples, file.path(out_dir, "samples_urine_only.csv"))
  message("Saved urine-only subset: ", file.path(out_dir, "samples_urine_only.csv"))
} else {
  message("No urine samples found.")
}

# ---- Plots for urine-only subset ----
if (nrow(urine_samples) > 0) {
  si_labels <- scales::label_number(scale_cut = scales::cut_si(" "))
  
  # 1) Top organisms (urine subset)
  if ("organism" %in% names(urine_samples)) {
    tab_org_u <- urine_samples %>%
      filter(!is.na(organism), organism != "") %>%
      count(organism, sort = TRUE) %>%
      slice_head(n = 20)
    
    if (nrow(tab_org_u) > 0) {
      p_org_u <- tab_org_u %>%
        mutate(organism = fct_reorder(organism, n)) %>%
        ggplot(aes(organism, n)) +
        geom_col() +
        coord_flip() +
        labs(title = "Urine-only subset: Top organisms",
             x = "Organism", y = "Count") +
        scale_y_continuous(labels = si_labels) +
        theme_minimal(base_size = 12)
      ggsave(file.path(out_dir, "UO1_counts_organism_urine_only.png"),
             p_org_u, width = 8, height = 6, dpi = 300)
      print(p_org_u)
    }
  }
  
  # 2) Host counts (urine subset)
  if ("host" %in% names(urine_samples)) {
    tab_host_u <- urine_samples %>%
      filter(!is.na(host), host != "") %>%
      count(host, sort = TRUE) %>%
      slice_head(n = 20)
    
    if (nrow(tab_host_u) > 0) {
      p_host_u <- tab_host_u %>%
        mutate(host = fct_reorder(host, n)) %>%
        ggplot(aes(host, n)) +
        geom_col() +
        coord_flip() +
        labs(title = "Urine-only subset: Hosts",
             x = "Host", y = "Count") +
        scale_y_continuous(labels = si_labels) +
        theme_minimal(base_size = 12)
      ggsave(file.path(out_dir, "UO2_counts_host_urine_only.png"),
             p_host_u, width = 8, height = 6, dpi = 300)
      print(p_host_u)
    }
  }
  
  # 3) Sample-source variants within urine subset (use the textual column you normalized)
  if ("sample_source" %in% names(urine_samples)) {
    tab_src_u <- urine_samples %>%
      filter(!is.na(sample_source), sample_source != "") %>%
      count(sample_source, sort = TRUE) %>%
      slice_head(n = 20)
    
    if (nrow(tab_src_u) > 0) {
      p_src_u <- tab_src_u %>%
        mutate(sample_source = fct_reorder(sample_source, n)) %>%
        ggplot(aes(sample_source, n)) +
        geom_col() +
        coord_flip() +
        labs(title = "Urine-only subset: Sample source labels",
             x = "Sample source", y = "Count") +
        scale_y_continuous(labels = si_labels) +
        theme_minimal(base_size = 12)
      ggsave(file.path(out_dir, "UO3_counts_sample_source_urine_only.png"),
             p_src_u, width = 8, height = 6, dpi = 300)
      print(p_src_u)
    }
  }
  
  # 4) Collection year histogram (urine subset)
  if ("collection_year" %in% names(urine_samples)) {
    p_year_uo <- urine_samples %>%
      filter(!is.na(collection_year), collection_year > 1900, collection_year < 2100) %>%
      ggplot(aes(x = collection_year)) +
      geom_histogram(binwidth = 1, boundary = 0, closed = "right") +
      scale_x_continuous(breaks = scales::pretty_breaks()) +
      labs(title = "Urine-only subset: Collection year",
           x = "Year", y = "Samples") +
      theme_minimal(base_size = 12)
    
    if (length(ggplot_build(p_year_uo)$data) > 0 && nrow(ggplot_build(p_year_uo)$data[[1]]) > 0) {
      ggsave(file.path(out_dir, "UO4_hist_year_urine_only.png"),
             p_year_uo, width = 8, height = 5, dpi = 300)
      print(p_year_uo)
    }
  }
  
  # 5) Stacked bar: Organism × Host (urine subset)
  if (all(c("organism","host") %in% names(urine_samples))) {
    tab_oh_uo <- urine_samples %>%
      filter(!is.na(organism), organism != "", !is.na(host), host != "") %>%
      count(organism, host, sort = TRUE)
    
    if (nrow(tab_oh_uo) > 0) {
      top_orgs_uo <- tab_oh_uo %>%
        group_by(organism) %>% summarise(n = sum(n), .groups = "drop") %>%
        slice_max(n, n = 15) %>% pull(organism)
      
      p_oh_uo <- tab_oh_uo %>%
        mutate(organism = if_else(organism %in% top_orgs_uo, organism, "Other")) %>%
        group_by(organism, host) %>% summarise(n = sum(n), .groups = "drop") %>%
        mutate(organism = fct_reorder(organism, n, sum)) %>%
        ggplot(aes(x = organism, y = n, fill = host)) +
        geom_col(position = "stack") +
        coord_flip() +
        labs(title = "Urine-only subset: Organism × Host",
             x = "Organism", y = "Count", fill = "Host") +
        theme_minimal(base_size = 12)
      ggsave(file.path(out_dir, "UO5_stack_organism_by_host_urine_only.png"),
             p_oh_uo, width = 9, height = 7, dpi = 300)
      print(p_oh_uo)
    }
  }
  
  # 6) Location  — optional
  if ("location_raw" %in% names(urine_samples)) {
    country_guess_uo <- urine_samples %>%
      mutate(country = location_raw) %>%
      mutate(country = str_replace_all(country, "\\s+", " ")) %>%
      mutate(country = str_replace(country, ".*:\\s*", "")) %>%
      mutate(country = str_split(country, ";|,") %>% purrr::map_chr(~ stringr::str_squish((.x[length(.x)] %||% NA_character_)))) %>%
      mutate(country = if_else(country == "", NA_character_, country))
    
    tab_loc <- country_guess_uo %>%
      filter(!is.na(country), country != "") %>%
      count(country, sort = TRUE) %>%
      slice_head(n = 30)
    
    if (nrow(tab_loc) > 0) {
      p_loc_uo <- tab_loc %>%
        mutate(country = fct_reorder(country, n)) %>%
        ggplot(aes(country, n)) +
        geom_col() +
        coord_flip() +
        labs(title = "Urine-only subset: Samples by location",
             x = "Location / Country", y = "Count") +
        theme_minimal(base_size = 12)
      ggsave(file.path(out_dir, "UO6_counts_location_urine_only.png"),
             p_loc_uo, width = 8, height = 6, dpi = 300)
      print(p_loc_uo)
    }
  }
} else {
  message("No urine samples to plot.")
}
