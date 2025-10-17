#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) {
    install.packages("optparse", repos = "https://cloud.r-project.org")
  }
  if (!requireNamespace("tidyverse", quietly = TRUE)) {
    install.packages("tidyverse", repos = "https://cloud.r-project.org")
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    install.packages("scales", repos = "https://cloud.r-project.org")
  }
  library(optparse)
  library(tidyverse)   # dplyr, readr, ggplot2, stringr, forcats, tibble, purrr
  library(scales)
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
              help="Comma-separated source terms, e.g. 'urine,blood' [default: %default]"),
  make_option(c("-H","--host"),   type="character", default="human",
              help="Comma-separated host filter, e.g. 'human' or 'human,mouse' [default: %default]"),
  make_option(c("-C","--search-extra-cols"), type="character", default="",
              help="Comma-separated columns to ALSO search (e.g., 'title,notes').")
)
opt <- parse_args(OptionParser(option_list = option_list))

in_csv  <- opt$input
out_dir <- opt$outdir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Parse source terms
source_terms <- stringr::str_split(opt$source, ",", simplify = TRUE) |> as.character() |> trimws()
source_terms <- source_terms[source_terms != ""]
if (length(source_terms) == 0) source_terms <- c("urine")
canon_terms <- unique(tolower(source_terms))

# Extra columns to search (optional)
extra_cols <- stringr::str_split(opt$`search-extra-cols`, ",", simplify = TRUE) |> as.character() |> trimws()
extra_cols <- extra_cols[extra_cols != ""]

# Hosts
host_terms <- stringr::str_split(opt$host, ",", simplify = TRUE) |> as.character() |> trimws()
host_terms <- host_terms[host_terms != ""]
if (length(host_terms) == 0) host_terms <- c("human")

# Scales helper
si_labels <- scales::label_number(scale_cut = scales::cut_short_scale())

# ---------------------------
# Helpers
# ---------------------------
pick_col <- function(df, candidates) {
  nms <- names(df)
  for (pat in candidates) {
    idx <- which(stringr::str_detect(nms, regex(pat, ignore_case = TRUE)))
    if (length(idx) > 0) return(nms[idx[1]])
  }
  return(NA_character_)
}

standardize_biosample <- function(df) {
  get_col <- function(patterns) {
    col <- pick_col(df, patterns)
    if (is.na(col)) return(rep(NA, nrow(df)))
    df[[col]]
  }
  tibble::tibble(
    accession        = get_col(c("^accession$","^acc$","^biosample[_ ]?id$","^sample[_ ]?accession$")),
    organism         = get_col(c("^organism$","organism[_ .]?name","^taxon","^species$","scientific[_ .]?name")),
    host             = get_col(c("^host$","host.*name","host.*species","host.*scientific.*name")),
    sample_source    = get_col(c("^sample[_ ]?source$","^isolation[_ ]?source$","^source$","specimen.*type","^material$","body[_ ]?site","anatomical")),
    collection_year  = suppressWarnings(as.integer(get_col(c("^collection[_ ]?year$","collection.*date","year"))))
  )
}

# Robust organism label (works across many header variants)
build_organism_label <- function(df) {
  candidates <- c("^organism$", "organism[_ .]?name", "scientific[_ .]?name", "^species$",
                  "taxon[_ .]?name", "tax[_ .]?name", "^taxon$", "^tax_id$", "taxonomy")
  nms <- names(df)
  for (pat in candidates) {
    idx <- which(stringr::str_detect(nms, regex(pat, ignore_case = TRUE)))
    if (length(idx) > 0) {
      v <- df[[ nms[idx[1]] ]]
      vch <- stringr::str_squish(as.character(v))
      if (any(!is.na(vch) & nzchar(vch))) return(vch)
    }
  }
  rep(NA_character_, nrow(df))
}

# Collapse to species: normalize strains/serovars/ST/etc. → Genus species
canon_species <- function(x) {
  if (length(x) == 0) return(character())
  y <- tolower(as.character(x))
  y <- stringr::str_replace_all(y, "_", " ")
  y <- stringr::str_squish(y)
  
  # Expand common genus abbreviations
  abbrev <- c(
    "^e\\.?\\s*coli\\b"               = "escherichia coli",
    "^k\\.?\\s*pneumoniae\\b"         = "klebsiella pneumoniae",
    "^p\\.?\\s*aeruginosa\\b"         = "pseudomonas aeruginosa",
    "^a\\.?\\s*baumannii\\b"          = "acinetobacter baumannii",
    "^s\\.?\\s*aureus\\b"             = "staphylococcus aureus",
    "^s\\.?\\s*pneumoniae\\b"         = "streptococcus pneumoniae",
    "^s\\.?\\s*pyogenes\\b"           = "streptococcus pyogenes",
    "^s\\.?\\s*enterica\\b"           = "salmonella enterica",
    "^e\\.?\\s*faecium\\b"            = "enterococcus faecium",
    "^e\\.?\\s*faecalis\\b"           = "enterococcus faecalis",
    "^h\\.?\\s*influenzae\\b"         = "haemophilus influenzae",
    "^n\\.?\\s*meningitidis\\b"       = "neisseria meningitidis",
    "^v\\.?\\s*cholerae\\b"           = "vibrio cholerae",
    "^l\\.?\\s*monocytogenes\\b"      = "listeria monocytogenes",
    "^m\\.?\\s*tuberculosis\\b"       = "mycobacterium tuberculosis",
    "^c\\.?\\s*difficile\\b"          = "clostridioides difficile"
  )
  for (pat in names(abbrev)) y <- stringr::str_replace(y, regex(pat), abbrev[[pat]])
  
  # Strip strain/serovar/ST descriptors and anything after
  y <- stringr::str_replace(y, "\\b(strain|substrain|isolate|clone|sequence type|sequence-type|st|clade|serovar|serotype|biovar|phylotype)\\b.*$", "")
  # Remove antigen codes like O157:H7, etc.
  y <- stringr::str_replace_all(y, "\\b(o|h|k)[0-9]+(?::[a-z0-9]+)?\\b", "")
  # Remove parentheses/brackets content
  y <- stringr::str_replace_all(y, "\\([^)]*\\)", "")
  y <- stringr::str_squish(y)
  
  # Normalize "Genus sp."
  y <- stringr::str_replace(y, "^(\\w+)\\s+sp\\.?\\b.*$", "\\1 sp.")
  
  toks <- stringr::str_split(y, "\\s+", simplify = TRUE)
  out <- character(length(y))
  for (i in seq_along(y)) {
    tks <- toks[i, ]; tks <- tks[tks != ""]
    if (length(tks) >= 2 && tks[2] != "sp.") {
      out[i] <- paste(tks[1], tks[2])
    } else if (length(tks) >= 1 && nzchar(tks[1])) {
      out[i] <- paste(tks[1], "sp.")
    } else out[i] <- NA_character_
  }
  out <- ifelse(is.na(out), NA_character_,
                paste0(stringr::str_to_title(word(out, 1)),
                       " ",
                       ifelse(word(out, 2) == "Sp.", "sp.", tolower(word(out, 2)))))
  out
}

# ---------------------------
# Read & prep (Option A: keep guessing, log parse problems)
# ---------------------------
message("Reading: ", in_csv)
raw <- suppressMessages(
  readr::read_csv(
    in_csv,
    show_col_types = FALSE,
    progress = FALSE,
    guess_max = 100000,
    locale = readr::locale(encoding = "UTF-8")
  )
)
pb <- readr::problems(raw)
if (inherits(pb, "tbl_df") && nrow(pb) > 0) {
  out_pb <- file.path(out_dir, "_read_problems.csv")
  readr::write_csv(pb, out_pb)
  message("Parsing issues detected (", nrow(pb), " rows). Details: ", out_pb)
} else {
  message("No parsing problems reported by readr.")
}

df_exp <- raw
std <- standardize_biosample(raw)
for (col in names(std)) if (!col %in% names(df_exp)) df_exp[[col]] <- std[[col]]

# Isolation explicit column
col_iso <- pick_col(raw, c("^isolation[_ ]?source$","isolation source","sample_source_detail",
                           "specimen_source","sample_type","specimen","material",
                           "body[_ ]?site","anatomical"))
df_exp <- df_exp %>%
  mutate(
    isolation_source_explicit = if (!is.na(col_iso)) raw[[col_iso]] else sample_source,
    isolation_source_explicit = stringr::str_squish(as.character(isolation_source_explicit))
  )

# Organism labels & species-collapsed label
df_exp$organism_label   <- build_organism_label(df_exp)
df_exp$organism_species <- canon_species(df_exp$organism_label)

# ---------------------------
# Host filter (default human)
# ---------------------------
host_cols <- c("host","host_name","host_common_name","host_scientific_name","host_species",
               "host_taxon","host_taxid","isolation_host","organism_host","host_description")
host_cols <- intersect(host_cols, names(df_exp))
if (length(host_cols) == 0) {
  host_cols <- names(df_exp)[stringr::str_detect(names(df_exp), regex("host", ignore_case = TRUE))]
}

per_host_regex <- setNames(vector("character", length(host_terms)), tolower(host_terms))
for (i in seq_along(host_terms)) {
  ht <- tolower(host_terms[[i]])
  if (ht == "human") {
    per_host_regex[[ht]] <- "(?i)human|homo\\s*sapiens|h\\.?\\s*sapiens|homo_sapiens|patient"
  } else {
    esc <- stringr::str_replace_all(ht, "([\\^$.|()?*+\\[\\]{}\\\\])", "\\\\\\1")
    per_host_regex[[ht]] <- paste0("(?i)", esc)
  }
}
host_combined_regex <- paste(per_host_regex, collapse = "|")

if (length(host_cols) > 0) {
  df_exp <- df_exp %>%
    mutate(across(all_of(host_cols), ~ tolower(as.character(.x)))) %>%
    filter(if_any(all_of(host_cols),
                  ~ stringr::str_detect(coalesce(.x, ""), regex(host_combined_regex))))
  message("Host filter kept rows: ", nrow(df_exp), " (hosts matched: ", paste(host_terms, collapse=", "), ")")
} else {
  message("No host-like columns found; host filter disabled.")
}

# ---------------------------
# Build search columns (NO title by default; add via -C if desired)
# ---------------------------
search_cols <- c(
  "sample_source","isolation_source_explicit","isolation_source","sample_source_detail",
  "specimen_source","sample_type","specimen","material","body_site","anatomical"
)
search_cols <- intersect(search_cols, names(df_exp))
auto_cols <- names(df_exp)[stringr::str_detect(names(df_exp), regex("(source|specimen)", ignore_case = TRUE))]
search_cols <- unique(c(search_cols, auto_cols))
if (length(extra_cols) > 0) {
  added <- intersect(extra_cols, names(df_exp))
  if (length(added) > 0) {
    search_cols <- unique(c(search_cols, added))
    message("Added to search columns: ", paste(added, collapse = ", "))
  } else {
    message("No valid extra search columns found among: ", paste(extra_cols, collapse = ", "))
  }
}
if (length(search_cols) == 0) {
  search_cols <- names(df_exp)[vapply(df_exp, function(v) is.character(v) || is.factor(v), logical(1))]
  message("Falling back to all character/factor columns as search space (", length(search_cols), ").")
}

# ---------------------------
# Regex per term
# ---------------------------
per_term_regex <- setNames(vector("character", length(canon_terms)), canon_terms)
for (i in seq_along(canon_terms)) {
  tm <- canon_terms[[i]]
  if (tm == "blood") {
    per_term_regex[[tm]] <- "(?i)blood"
  } else if (tm == "urine") {
    per_term_regex[[tm]] <- "(?i)(urine|urinary|uti)"
  } else if (tm %in% c("feces","faeces","stool","fecal","faecal")) {
    per_term_regex[[tm]] <- "(?i)(feces|faeces|stool|fecal|faecal|rectal|perianal|peri-anal)"
  } else {
    esc <- stringr::str_replace_all(tm, "([\\^$.|()?*+\\[\\]{}\\\\])", "\\\\\\1")
    per_term_regex[[tm]] <- paste0("(?i)", esc)
  }
}
combined_regex <- paste(per_term_regex, collapse = "|")

# ---------------------------
# Subsets
# ---------------------------
subsets_by_term <- list()
term_hits_summaries <- list()

# Combined subset
combined_subset <- df_exp %>%
  mutate(across(all_of(search_cols), ~ tolower(as.character(.x)))) %>%
  filter(if_any(all_of(search_cols), ~ str_detect(coalesce(.x, ""), regex(combined_regex)))) %>%
  distinct(accession, .keep_all = TRUE)

subset_tag <- tolower(paste(canon_terms, collapse = "_")) |> gsub("\\W+","_", x = _)
if (nrow(combined_subset) > 0) {
  out_comb <- file.path(out_dir, paste0("samples_", subset_tag, "_only.csv"))
  readr::write_csv(combined_subset, out_comb)
  message("Saved combined subset (", paste(canon_terms, collapse=","), "): ", out_comb)
} else {
  message("No rows matched any of: ", paste(canon_terms, collapse=", "))
}

# Per-term subsets + plots
for (term in canon_terms) {
  term_regex <- per_term_regex[[term]]
  subset_term <- df_exp %>%
    mutate(across(all_of(search_cols), ~ tolower(as.character(.x)))) %>%
    filter(if_any(all_of(search_cols), ~ str_detect(coalesce(.x, ""), regex(term_regex)))) %>%
    distinct(accession, .keep_all = TRUE)
  
  term_tag <- tolower(gsub("\\W+","_", term))
  subsets_by_term[[term_tag]] <- subset_term
  message("Subset size for '", term, "': ", nrow(subset_term))
  
  # Diagnostics per column
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
    
    # --- Per-term plots (species-collapsed) ---
    if ("organism_species" %in% names(subset_term)) {
      tab_sp <- subset_term %>%
        filter(!is.na(organism_species), organism_species != "") %>%
        mutate(organism_species = as.character(organism_species)) %>%
        count(organism_species, sort = TRUE) %>%
        slice_head(n = 20)
      if (nrow(tab_sp) > 0) {
        p1 <- tab_sp %>%
          mutate(organism_species = forcats::fct_reorder(organism_species, n)) %>%
          ggplot(aes(organism_species, n)) +
          geom_col() + coord_flip() +
          labs(title = paste("Subset:", term, "— Top species (strain-collapsed)"),
               x = "Species", y = "Count") +
          scale_y_continuous(labels = si_labels) +
          theme_minimal(base_size = 12)
        ggsave(file.path(out_dir, paste0("SUB1_counts_species_", term_tag, "_only.png")),
               p1, width = 8, height = 6, dpi = 300)
      }
    }
    
    if ("host" %in% names(subset_term)) {
      tab_host <- subset_term %>%
        filter(!is.na(host), host != "") %>%
        mutate(host = as.character(host)) %>%
        count(host, sort = TRUE) %>% slice_head(n = 20)
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
    
    if ("sample_source" %in% names(subset_term)) {
      tab_src <- subset_term %>%
        filter(!is.na(sample_source), sample_source != "") %>%
        mutate(sample_source = as.character(sample_source)) %>%
        count(sample_source, sort = TRUE) %>% slice_head(n = 20)
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
    
    # --- Species counts per year, per term (stacked) ---
    if (all(c("organism_species","collection_year") %in% names(subset_term))) {
      sp_year <- subset_term %>%
        filter(!is.na(collection_year), collection_year > 1900, collection_year < 2100,
               !is.na(organism_species), organism_species != "") %>%
        mutate(organism_species = as.character(organism_species)) %>%
        count(collection_year, organism_species, name = "n", sort = TRUE)
      
      if (nrow(sp_year) > 0) {
        # top species within this term across all years
        top_sp <- sp_year %>%
          group_by(organism_species) %>% summarise(N = sum(n), .groups = "drop") %>%
          slice_max(N, n = 10) %>% pull(organism_species)
        sp_year2 <- sp_year %>%
          mutate(organism_species = if_else(organism_species %in% top_sp, organism_species, "Other")) %>%
          group_by(collection_year, organism_species) %>% summarise(n = sum(n), .groups = "drop")
        
        p5 <- sp_year2 %>%
          mutate(organism_species = forcats::fct_reorder(organism_species, n, .fun = sum)) %>%
          ggplot(aes(x = collection_year, y = n, fill = organism_species)) +
          geom_col() +
          scale_x_continuous(breaks = scales::pretty_breaks()) +
          labs(title = paste("Subset:", term, "— Species counts per year (top 10 + Other)"),
               x = "Year", y = "Samples", fill = "Species") +
          theme_minimal(base_size = 12)
        ggsave(file.path(out_dir, paste0("SUB5_species_by_year_", term_tag, "_only.png")),
               p5, width = 10, height = 6, dpi = 300)
      }
    }
    
  } else {
    message("No samples matched source term: ", term)
  }
}

# ---------------------------
# COMBO comparison plots (any N terms >= 2)
# ---------------------------
TOP_K_GLOBAL   <- 15   # shared top-K species globally
TOP_K_PER_TERM <- 10   # per-term top-K
YEAR_BINWIDTH  <- 1    # years per bin

avail_terms <- names(subsets_by_term)
have_any <- avail_terms[vapply(subsets_by_term[avail_terms], function(df) !is.null(df) && nrow(df) > 0, logical(1))]

if (length(have_any) >= 2) {
  combined_df_all <- dplyr::bind_rows(lapply(have_any, function(tg) {
    df <- subsets_by_term[[tg]]
    if (!is.null(df) && nrow(df) > 0) { df$source_term <- tg; df } else { NULL }
  }))
  
  if (!is.null(combined_df_all) && nrow(combined_df_all) > 0) {
    term_sizes <- combined_df_all %>% dplyr::count(source_term, name = "N")
    
    # SPECIES (strain-collapsed): counts + proportions + per-term top-K
    if (all(c("organism_species","source_term") %in% names(combined_df_all))) {
      sp_raw <- combined_df_all %>%
        filter(!is.na(organism_species), organism_species != "") %>%
        mutate(organism_species = as.character(organism_species))
      
      if (nrow(sp_raw) > 0) {
        sp_counts <- sp_raw %>% count(source_term, organism_species, name = "n", sort = TRUE)
        
        top_sp_all <- sp_counts %>%
          group_by(organism_species) %>% summarise(N = sum(n), .groups = "drop") %>%
          arrange(desc(N)) %>% slice_head(n = TOP_K_GLOBAL) %>% pull(organism_species) %>% as.character()
        
        sp_shared <- sp_counts %>%
          mutate(organism_species = ifelse(organism_species %in% top_sp_all, organism_species, "Other")) %>%
          group_by(source_term, organism_species) %>% summarise(n = sum(n), .groups = "drop")
        
        if (nrow(sp_shared) > 0) {
          # counts
          p_sp_counts <- sp_shared %>%
            mutate(organism_species = ifelse(is.na(organism_species) | organism_species == "", "Unknown", organism_species)) %>%
            mutate(organism_species = forcats::fct_reorder(organism_species, n, .fun = sum)) %>%
            ggplot(aes(x = organism_species, y = n, fill = source_term)) +
            geom_col(position = "dodge") + coord_flip() +
            labs(title = "Species by source term (strain-collapsed, shared top-K)",
                 x = "Species", y = "Count", fill = "Source term") +
            scale_y_continuous(labels = si_labels) + theme_minimal(base_size = 12)
          ggsave(file.path(out_dir, "COMBO_counts_species_all_terms.png"),
                 p_sp_counts, width = 10, height = 8, dpi = 300)
          
          # proportions
          sp_props <- sp_shared %>% group_by(source_term) %>% mutate(prop = n / sum(n)) %>% ungroup()
          if (nrow(sp_props) > 0) {
            p_sp_props <- sp_props %>%
              mutate(organism_species = forcats::fct_reorder(organism_species, prop, .fun = sum)) %>%
              ggplot(aes(x = organism_species, y = prop, fill = source_term)) +
              geom_col(position = "dodge") + coord_flip() +
              scale_y_continuous(labels = label_percent(accuracy = 1)) +
              labs(title = "Species by source term (proportions)",
                   x = "Species", y = "Share", fill = "Source term") +
              theme_minimal(base_size = 12)
            ggsave(file.path(out_dir, "COMBO_props_species_all_terms.png"),
                   p_sp_props, width = 10, height = 8, dpi = 300)
          }
          
          # per-term top-K
          sp_perterm <- sp_counts %>% group_by(source_term) %>%
            slice_max(n, n = TOP_K_PER_TERM, with_ties = FALSE) %>% ungroup()
          if (nrow(sp_perterm) > 0) {
            p_sp_perterm <- sp_perterm %>%
              mutate(organism_species = as.character(organism_species)) %>%
              group_by(source_term) %>%
              mutate(organism_species = forcats::fct_reorder(organism_species, n)) %>%
              ggplot(aes(x = organism_species, y = n, fill = source_term)) +
              geom_col(show.legend = FALSE) + coord_flip() +
              facet_wrap(~ source_term, scales = "free_y") +
              labs(title = "Species (top-K within each term)",
                   x = "Species", y = "Count") +
              scale_y_continuous(labels = si_labels) +
              theme_minimal(base_size = 12)
            ggsave(file.path(out_dir, "COMBO_counts_species_per_term.png"),
                   p_sp_perterm, width = 12, height = 8, dpi = 300)
          }
        } else {
          message("COMBO (species): no rows after top-K bucketing; skipped.")
        }
      } else {
        message("COMBO (species): no organism_species values; skipped.")
      }
    } else {
      message("COMBO (species): required columns missing; skipped.")
    }
    
    # YEARS (counts + proportions; with Ns)
    if (all(c("collection_year","source_term") %in% names(combined_df_all))) {
      year_df <- combined_df_all %>% filter(!is.na(collection_year), collection_year > 1900, collection_year < 2100)
      
      # counts
      p_year_counts <- year_df %>%
        ggplot(aes(x = collection_year)) +
        geom_histogram(binwidth = YEAR_BINWIDTH, boundary = 0, closed = "right") +
        facet_wrap(~ source_term, ncol = 1, scales = "free_y",
                   labeller = labeller(source_term = function(x) {
                     paste0(x, " (N=", (term_sizes$N[match(x, term_sizes$source_term)]), ")")
                   })) +
        scale_x_continuous(breaks = pretty_breaks()) +
        labs(title = "Collection year by source term (counts)",
             subtitle = paste0("Binwidth = ", YEAR_BINWIDTH, " year"),
             x = "Year", y = "Samples") +
        theme_minimal(base_size = 12)
      ggsave(file.path(out_dir, "COMBO_hist_year_all_terms.png"),
             p_year_counts, width = 10, height = 10, dpi = 300)
      
      # proportions
      p_year_props <- year_df %>%
        group_by(source_term, collection_year) %>% summarise(n = n(), .groups = "drop") %>%
        group_by(source_term) %>% mutate(prop = n / sum(n)) %>% ungroup() %>%
        ggplot(aes(x = collection_year, y = prop)) +
        geom_col() +
        facet_wrap(~ source_term, ncol = 1,
                   labeller = labeller(source_term = function(x) {
                     paste0(x, " (N=", (term_sizes$N[match(x, term_sizes$source_term)]), ")")
                   })) +
        scale_y_continuous(labels = label_percent(accuracy = 1)) +
        scale_x_continuous(breaks = pretty_breaks()) +
        labs(title = "Collection year by source term (proportions)",
             x = "Year", y = "Share of term") +
        theme_minimal(base_size = 12)
      ggsave(file.path(out_dir, "COMBO_hist_year_props_all_terms.png"),
             p_year_props, width = 10, height = 10, dpi = 300)
    } else {
      message("COMBO (year): required columns missing; skipped.")
    }
    
    # --- COMBO species counts per year (stacked, faceted by term) ---
    if (all(c("organism_species","collection_year","source_term") %in% names(combined_df_all))) {
      sp_year_all <- combined_df_all %>%
        filter(!is.na(collection_year), collection_year > 1900, collection_year < 2100,
               !is.na(organism_species), organism_species != "") %>%
        mutate(organism_species = as.character(organism_species)) %>%
        count(source_term, collection_year, organism_species, name = "n", sort = TRUE)
      
      if (nrow(sp_year_all) > 0) {
        # top species per term (across years)
        top_per_term <- sp_year_all %>%
          group_by(source_term, organism_species) %>% summarise(N = sum(n), .groups = "drop") %>%
          group_by(source_term) %>% slice_max(N, n = 10, with_ties = FALSE)
        
        sp_year2 <- sp_year_all %>%
          left_join(top_per_term %>% mutate(top_flag = TRUE),
                    by = c("source_term","organism_species")) %>%
          mutate(organism_species = if_else(!is.na(top_flag) & top_flag, organism_species, "Other")) %>%
          group_by(source_term, collection_year, organism_species) %>%
          summarise(n = sum(n), .groups = "drop")
        
        p_sy_combo <- sp_year2 %>%
          group_by(source_term) %>%
          mutate(organism_species = forcats::fct_reorder(organism_species, n, .fun = sum)) %>%
          ungroup() %>%
          ggplot(aes(x = collection_year, y = n, fill = organism_species)) +
          geom_col() +
          facet_wrap(~ source_term, ncol = 1, scales = "free_y") +
          scale_x_continuous(breaks = scales::pretty_breaks()) +
          labs(title = "Species counts per year by source term (top 10 per term + Other)",
               x = "Year", y = "Samples", fill = "Species") +
          theme_minimal(base_size = 12)
        ggsave(file.path(out_dir, "COMBO_species_by_year_all_terms.png"),
               p_sy_combo, width = 11, height = 10, dpi = 300)
      } else {
        message("COMBO (species-by-year): nothing to plot; skipped.")
      }
    } else {
      message("COMBO (species-by-year): required columns missing; skipped.")
    }
    
  } else {
    message("COMBO: combined_df_all had 0 rows; skipped.")
  }
} else {
  message("COMBO: need at least 2 non-empty term subsets; skipped.")
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

# ---------------------------
# Final summary of CSV outputs (no parse warnings)
# ---------------------------
message("\nFinal CSV table dimensions:")

# Build a small registry from the objects we already have in memory
dim_lines <- list()

# Combined
if (exists("combined_subset") && is.data.frame(combined_subset)) {
  dim_lines <- append(dim_lines, sprintf(
    "samples_%s_only.csv: %d rows × %d cols",
    subset_tag, nrow(combined_subset), ncol(combined_subset)
  ))
}

# Per-term
for (term in canon_terms) {
  term_tag <- tolower(gsub("\\W+","_", term))
  df <- subsets_by_term[[term_tag]]
  if (is.data.frame(df) && nrow(df) > 0) {
    dim_lines <- append(dim_lines, sprintf(
      "samples_%s_only.csv: %d rows × %d cols",
      term_tag, nrow(df), ncol(df)
    ))
  }
}

# If for any reason a file exists but we didn't have the df (or it had 0 rows),
# list it and read as all-character to avoid parse warnings, just for dimensions.
all_csvs <- list.files(out_dir, pattern = "^samples_.*\\.csv$", full.names = TRUE)
reported <- sub("^samples_(.*)\\.csv$", "samples_\\1.csv", basename(gsub("_only\\.csv$","_only.csv", gsub(".*?/","", sub("^.*/","", all_csvs)))))
already <- sub("^samples_(.*)$", "samples_\\1", gsub("\\.csv$","", sub("^.*?/","", dim_lines)))

still_need <- setdiff(basename(all_csvs), gsub("^.*?:\\s+","", paste0(sub("^(samples_[^:]+):.*$","\\1", dim_lines), collapse=";")))
if (length(still_need) > 0) {
  for (f in all_csvs) {
    if (basename(f) %in% still_need) {
      df <- tryCatch(
        suppressWarnings(readr::read_csv(
          f,
          col_types = readr::cols(.default = readr::col_character()),
          show_col_types = FALSE,
          progress = FALSE
        )),
        error = function(e) NULL
      )
      if (!is.null(df)) {
        dim_lines <- append(dim_lines, sprintf(
          "%s: %d rows × %d cols", basename(f), nrow(df), ncol(df)
        ))
      } else {
        dim_lines <- append(dim_lines, sprintf("%s: [could not read]", basename(f)))
      }
    }
  }
}

# Print in a stable order (combined first, then alphabetical)
dim_lines <- unique(dim_lines)
dim_lines <- dim_lines[order(tolower(dim_lines))]

purrr::walk(dim_lines, message)  