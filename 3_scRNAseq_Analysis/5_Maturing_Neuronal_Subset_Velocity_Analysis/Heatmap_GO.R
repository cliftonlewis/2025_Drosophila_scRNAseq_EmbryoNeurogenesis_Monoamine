# ------------------------ setup ------------------------
options(stringsAsFactors = FALSE)
set.seed(42)

# Adjust to your project root
setwd("/DataDrives/Drive2/Clifton/R_Projects/2025_Drosophila_scRNAseq_MonoamineSpecification/ANALYSIS/Step_3_scVelo/Step_3.2_Mature_neuronal_subset_analysis/")

# ---- Packages ----
library(tidyverse)
library(scales)

# ------------------------ config ------------------------
input_dir  <- "/DataDrives/Drive2/Clifton/R_Projects/2025_Drosophila_scRNAseq_MonoamineSpecification/ANALYSIS/Step_3_scVelo/Step_3.2_Mature_neuronal_subset_analysis/"
file_regex <- "^GO_bin_\\d+\\.csv$"  # e.g., GO_bin_0.csv ... GO_bin_19.csv
N          <- 5                      # top terms per bin to retain in the union
bins_all   <- 0:19                   # expected latent-time bins
cap        <- 10                     # cap for -log10(p) in the color scale
wrap_width <- 40

# ------------------------ files ------------------------
files <- list.files(input_dir, pattern = file_regex, full.names = TRUE)
stopifnot(length(files) > 0)

# ------------------------ robust reader ------------------------
read_one_bin <- function(f) {
  bin <- as.integer(stringr::str_match(basename(f), "GO_bin_(\\d+)\\.csv$")[,2])
  df  <- readr::read_csv(f, show_col_types = FALSE)
  
  # g:Profiler columns vary; find likely term & p columns
  term_col <- intersect(c("name","term_name","term"), names(df))[1]
  p_col    <- intersect(c("p_value","p.value","adjusted_p_value","padj"), names(df))[1]
  if (is.na(term_col) || is.na(p_col)) stop("Missing term/p columns in ", f)
  
  # Parse p robustly: handles "1e-6", "0.001", "<1e-300" (becomes 1e-300), etc.
  p_raw <- as.character(df[[p_col]])
  p_num <- suppressWarnings(readr::parse_number(p_raw, locale = readr::locale(decimal_mark = ".", grouping_mark = ",")))
  
  # If there were strings like "<1e-300>", set those to a tiny floor
  is_lt <- grepl("^\\s*<\\s*", p_raw)
  p_num[is_lt & is.na(p_num)] <- 1e-300
  
  # If any zeros slipped in, keep as 0; -log10(0) -> Inf, later capped to 'cap'
  tibble(
    bin   = bin,
    term  = as.character(df[[term_col]]),
    p_val = as.numeric(p_num),
    score = -log10(p_num)
  )
}

go <- purrr::map_dfr(files, read_one_bin)

# ------------------------ rank & keep union of top-N ------------------------
go_ranked <- go %>%
  group_by(bin) %>%
  arrange(desc(score), .by_group = TRUE) %>%
  mutate(rank_in_bin = row_number(),
         topN = rank_in_bin <= N) %>%
  ungroup()

keep_terms <- go_ranked %>% filter(topN) %>% distinct(term)

# Keep all bins for these terms (fill missing with NA), and keep topN flags per cell
go_keep <- go %>%
  semi_join(keep_terms, by = "term") %>%
  group_by(term) %>%
  tidyr::complete(bin = bins_all, fill = list(p_val = NA_real_, score = NA_real_)) %>%
  ungroup() %>%
  left_join(go_ranked %>% select(term, bin, topN), by = c("term","bin")) %>%
  mutate(topN = tidyr::replace_na(topN, FALSE))

# ------------------------ row order: bin of peak enrichment ------------------------
row_order <- go_keep %>%
  group_by(term) %>%
  summarize(
    peak_bin = bin[which.max(replace_na(score, -Inf))][1],  # leftmost if tie
    .groups = "drop"
  ) %>%
  arrange(peak_bin, term)

# ------------------------ plot data ------------------------
go_plot <- go_keep %>%
  mutate(
    bin          = factor(bin, levels = bins_all),
    term         = factor(term, levels = rev(row_order$term)),  # reverse for top at top
    score_capped = pmin(score, cap),
    sig = case_when(
      is.na(p_val) ~ "",
      p_val < 1e-4 ~ "****",
      p_val < 1e-3 ~ "***",
      p_val < 1e-2 ~ "**",
      p_val < 5e-2 ~ "*",
      TRUE         ~ ""
    ),
    # factor used only to make a legend; empty strings -> NA to hide
    sig_legend = factor(
      dplyr::na_if(sig, ""),
      levels = c("****","***","**","*")
    )
  )

# ------------------------ plot ------------------------
p <- ggplot(go_plot, aes(x = bin, y = term, fill = score_capped)) +
  geom_tile() +
  
  # Invisible points to create a legend for significance thresholds
  geom_point(
    data = subset(go_plot, topN & !is.na(sig_legend)),
    aes(shape = sig_legend),
    size = 0, alpha = 0
  ) +
  
  # Asterisks on significant top-N tiles
  geom_text(
    data = subset(go_plot, topN & sig != ""),
    aes(label = sig),
    color = "white", size = 3, fontface = "bold"
  ) +
  
  scale_fill_gradientn(
    colours  = c("#FFE39A", "#F48FB1", "#7B1FA2"),   # yellow → pink → purple
    values   = scales::rescale(c(0, cap/2, cap)),
    na.value = "black",                               # missing shown as black
    limits   = c(0, cap),
    name     = expression(-log[10](p))
  ) +
  
  # Legend for significance thresholds (stars)
  scale_shape_manual(
    name   = "Significance (p)",
    values = c(8, 8, 8, 8),  # star symbol in legend
    labels = c("< 1e-4 (****)", "< 1e-3 (***)", "< 1e-2 (**)","< 0.05 (*)")
  ) +
  guides(
    shape = guide_legend(
      override.aes = list(size = 4, alpha = 1, colour = "black")
    )
  ) +
  
  scale_y_discrete(labels = \(x) stringr::str_wrap(x, width = wrap_width)) +
  labs(
    x = "Latent time bins (0–19)",
    y = sprintf("GO enrichment terms (union of top %d/bin)", N),
    title = "Neurogenic Lineage – Dynamical Genes – GO"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid  = element_blank(),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 9),
    plot.title  = element_text(face = "bold", size = 14)
  )

# Show in device
p

# ------------------------ save (optional) ------------------------
ggsave(file.path(input_dir, "GO_heatmap_termsY.png"), p, width = 7.5, height = 10, dpi = 300)
ggsave(file.path(input_dir, "GO_heatmap_by_latent_time_vertical.pdf"), p, width = 10, height = 8)
ggsave(file.path(input_dir, "GO_heatmap_by_latent_time_vertical.png"), p, width = 10, height = 8, dpi = 300)
