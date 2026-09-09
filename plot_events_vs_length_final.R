# ============================================================
# Events (raw count) vs. gene alignment length (raw bp) -- SNP,
# Frameshift, and In-frame indel, Core and Shell separately.
#
# Reverted from the events-per-kb rate version back to raw counts
# (rate version risked the small-denominator/ratio-correlation
# artifact flagged earlier -- raw count avoids that entirely).
#
# Styling polished to match this project's established figure
# conventions: no gridlines, L-shaped axes only, refined point
# styling, clean R² annotation, consistent typography.
# ============================================================

library(tidyverse)
library(ggpmisc)
library(ggrepel)

base_dir <- "~/Downloads/macrogen_1c/pangenome/pangenome_results"
category_palette <- c("Core" = "#B491E5", "Shell" = "#66A8E5")

theme_figure <- function(base_size = 12) {
  theme_minimal(base_size = base_size, base_family = "Helvetica") +
    theme(
      panel.grid          = element_blank(),
      axis.line.x         = element_line(color = "black", linewidth = 0.4),
      axis.line.y         = element_line(color = "black", linewidth = 0.4),
      axis.ticks           = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length    = unit(3, "pt"),
      axis.title            = element_text(face = "bold", size = base_size),
      axis.text              = element_text(size = base_size - 1, color = "black"),
      strip.background        = element_blank(),
      strip.text               = element_text(face = "bold", size = base_size + 0.5),
      legend.position            = "none",
      panel.spacing                = unit(1.6, "lines"),
      plot.title                     = element_text(face = "bold", size = base_size + 3),
      plot.subtitle                    = element_text(size = base_size - 1, color = "grey40"),
      plot.margin                        = margin(14, 18, 10, 10)
    )
}

# ── Base family list + misclustering exclusion ────────────────────
base_families <- read_csv(file.path(base_dir, "snp_indel_events.csv"), show_col_types = FALSE) %>%
  select(category, family_id, aln_len, snp)

screen <- read_csv(file.path(base_dir, "family_identity_screen.csv"), show_col_types = FALSE)
flagged <- screen %>% filter(flagged_misclustered == TRUE) %>%
  mutate(flag_key = paste(category, family_id)) %>% pull(flag_key)

base_families <- base_families %>%
  mutate(flag_key = paste(category, family_id)) %>%
  filter(!flag_key %in% flagged, category %in% c("Core", "Shell"))

large_effect <- read_csv(file.path(base_dir, "large_effect_detail.csv"), show_col_types = FALSE)

# ── Gene product names, for labeling outliers meaningfully ────────
syn_dir <- file.path(base_dir, "synteny")
genes <- read_csv(file.path(syn_dir, "genes.tsv"), show_col_types = FALSE)
products <- read_csv(file.path(syn_dir, "gene_products.csv"), show_col_types = FALSE)

family_products <- genes %>%
  filter(!is.na(family_id)) %>%
  left_join(products, by = "gene_id") %>%
  filter(!is.na(product)) %>%
  count(family_id, product, sort = TRUE) %>%
  group_by(family_id) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(family_id, product)

build_count_data <- function(event_label) {
  if (event_label == "SNP") {
    df <- base_families %>% mutate(n_events = snp)
  } else {
    event_type <- if (event_label == "Frameshift") "frameshift" else "in_frame_indel"
    counts <- large_effect %>% filter(event_type == !!event_type) %>%
      count(category, family_id, name = "n_events")
    df <- base_families %>%
      left_join(counts, by = c("category", "family_id")) %>%
      mutate(n_events = replace_na(n_events, 0))
  }
  df
}

make_plot <- function(data, event_label, out_filename, exclude_product_pattern = NULL) {
  rate_summary <- data %>%
    group_by(category) %>%
    summarise(mean_rate = sum(n_events) / sum(aln_len), n_families = n(), .groups = "drop")

  cat("\n===", event_label, "===\n")
  print(rate_summary)

  # ── Identify outliers by RESIDUAL from the fitted trend (per
  # category) -- the principled definition of "outlier" here: biggest
  # gap between actual and predicted count, in EITHER direction. Top 3
  # over-mutated + top 3 under-mutated per category, labeled with real
  # gene product names.
  data_with_resid <- data %>%
    group_by(category) %>%
    mutate(fitted = predict(lm(n_events ~ aln_len)), residual = n_events - fitted) %>%
    ungroup() %>%
    left_join(family_products, by = "family_id") %>%
    mutate(product = replace_na(product, "hypothetical protein"))

  outliers <- data_with_resid %>%
    group_by(category) %>%
    group_modify(~ bind_rows(
      slice_max(.x, residual, n = 3),
      slice_min(.x, residual, n = 3)
    )) %>%
    ungroup() %>%
    mutate(outlier_label = paste0(product, " (", family_id, ")"))

  if (!is.null(exclude_product_pattern)) {
    n_before <- nrow(outliers)
    outliers <- outliers %>% filter(!str_detect(product, regex(exclude_product_pattern, ignore_case = TRUE)))
    cat("Excluded", n_before - nrow(outliers), "outlier(s) matching '", exclude_product_pattern,
        "' -- label was ambiguous (unclear which point it pointed to)\n")
  }

  cat("\nOutliers (top 3 over- and under-mutated per category):\n")
  print(outliers %>% select(category, family_id, product, aln_len, n_events, residual))

  p <- ggplot(data, aes(x = aln_len, y = n_events)) +
    geom_point(aes(color = category), alpha = 0.4, size = 1.1, shape = 16) +
    geom_smooth(method = "lm", se = TRUE, color = "grey20", linewidth = 0.65,
                fill = "grey75", alpha = 0.25) +
    geom_abline(data = rate_summary, aes(slope = mean_rate, intercept = 0),
                linetype = "dashed", color = "grey45", linewidth = 0.55) +
    geom_point(data = outliers, shape = 21, size = 2.2, color = "black", fill = "white", stroke = 0.8) +
    geom_text_repel(data = outliers, aes(label = outlier_label),
                     size = 2.5, family = "Helvetica", fontface = "bold", color = "grey15",
                     max.overlaps = Inf, min.segment.length = 0, segment.color = "grey50",
                     segment.linewidth = 0.3, box.padding = 0.5, seed = 42) +
    stat_poly_eq(aes(label = paste(after_stat(rr.label))),
                 formula = y ~ x, size = 3.6, family = "Helvetica", fontface = "bold",
                 label.x = 0.04, label.y = 0.95, color = "grey20") +
    scale_color_manual(values = category_palette) +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::comma) +
    facet_wrap(~ category, scales = "free", ncol = 2) +
    labs(
      title = paste0(event_label, " events vs. gene length"),
      subtitle = "Dashed line = expected count under each category's mean rate  |  Labeled points = top 3 over/under-mutated outliers by regression residual",
      x = "Gene (alignment) length, bp",
      y = paste0(event_label, " events (raw count)")
    ) +
    theme_figure()

  out_path <- file.path(base_dir, out_filename)
  ggsave(out_path, p, width = 13, height = 6.5, dpi = 300, bg = "white")
  cat("Saved:", out_path, "\n")
}

snp_data        <- build_count_data("SNP")
frameshift_data <- build_count_data("Frameshift")
inframe_data    <- build_count_data("In-frame indel")

make_plot(snp_data, "SNP", "snp_vs_length_final.png", exclude_product_pattern = "methylmalonyl")
make_plot(frameshift_data, "Frameshift", "frameshift_vs_length_final.png")
make_plot(inframe_data, "In-frame indel", "inframe_indel_vs_length_final.png")
