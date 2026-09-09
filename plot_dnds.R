# ============================================================
# Evolutionary rate (dN, dS, dN/dS) per gene: Core vs Shell
# Three SEPARATE panels (independent y-scales -- dN, dS, and dN/dS
# are different quantities with different typical magnitudes)
#
# Same visual language as plot_snp_indel_events.R:
#   - Core = violet (Lavender/Purple), Shell = blue (Sky Blue).
#   - Bars: thin, flat solid fill, NO outline, tightly dodged at one
#     shared x-position per panel (not two separate x-axis categories).
#   - Data points: black-and-white (black outline, no fill); shape
#     distinguishes Core (circle) vs Shell (square). Minimum 10
#     representative points per group (nearest real data point to each
#     of the 10th/20th/.../100th percentiles, or all points if fewer
#     than 10 are available).
#   - Error bars: capped ("I" shape), black, floored at 0, drawn LAST
#     so they render frontmost (on top of the points), matching the
#     fix applied to plot_heaps_law.R.
#   - Legend re-enabled, no border/background.
#   - Significance shown as the actual p-value, plain line, no
#     drop-tick ends.
#   - Exclusions (misclustered families, non-codon-frame-safe
#     alignments) already happened upstream in calculate_dnds.py --
#     this script does not need to re-filter.
#   - dN/dS is undefined for families with mean dS = 0 (blank in the
#     source CSV); those rows are dropped specifically from the dN/dS
#     panel, not from dN or dS.
#   - Three independent statistical tests (dN, dS, dN/dS), each chosen
#     conditionally: Shapiro-Wilk normality check per group decides
#     Welch's t-test (both ~normal) vs. Wilcoxon rank-sum (otherwise).
# ============================================================

library(tidyverse)

# ── 0. Shared palette + theme ──────────────────────────────────
category_palette <- c("Core" = "#B491E5", "Shell" = "#66A8E5")  # violet, blue
category_shape   <- c("Core" = 1, "Shell" = 0)                   # open circle, open square
N_REPRESENTATIVE <- 10

theme_figure <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "Helvetica") +
    theme(
      panel.grid        = element_blank(),
      axis.line         = element_line(color = "black", linewidth = 0.4),
      axis.ticks        = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(2, "pt"),
      axis.title        = element_text(face = "bold", size = base_size),
      axis.text         = element_text(size = base_size - 1, color = "black"),
      strip.background   = element_blank(),
      strip.text         = element_text(face = "bold", size = base_size),
      panel.spacing      = unit(16, "pt"),
      legend.position     = "top",
      legend.title        = element_blank(),
      legend.text         = element_text(size = base_size - 1),
      legend.background   = element_blank(),
      legend.key          = element_blank()
    )
}

format_pval <- function(p) {
  if (p < 0.0001) return("p < 0.0001")
  paste0("p = ", formatC(p, format = "f", digits = 4))
}

select_representative_points <- function(data, n_points = N_REPRESENTATIVE) {
  data %>%
    group_by(metric, category) %>%
    group_modify(~ {
      d <- .x
      if (nrow(d) <= n_points) return(d)
      targets <- quantile(d$value, probs = seq(0, 1, length.out = n_points), type = 7)
      idx <- unique(sapply(targets, function(t) which.min(abs(d$value - t))))
      d[idx, ]
    }) %>%
    ungroup()
}

# ── 1. Load data (exclusions already applied upstream in Python) ──
base_dir <- "~/Downloads/macrogen_1c/pangenome/pangenome_results"
csv_file <- file.path(base_dir, "dnds_per_family.csv")

df_raw <- read_csv(csv_file, show_col_types = FALSE) %>%
  mutate(category = factor(category, levels = c("Core", "Shell")))

cat("Families per category (already excludes misclustered/frame-unsafe):\n")
print(table(df_raw$category))

df <- df_raw %>%
  select(category, family_id, dN = mean_dN, dS = mean_dS, `dN/dS` = dN_dS) %>%
  pivot_longer(c(dN, dS, `dN/dS`), names_to = "metric", values_to = "value") %>%
  filter(!is.na(value)) %>%   # drops undefined dN/dS (mean dS == 0) from that panel only
  mutate(metric = factor(metric, levels = c("dN", "dS", "dN/dS")))

cat("\nFamilies per metric x category (after dropping undefined dN/dS):\n")
print(table(df$metric, df$category))

# ── 2. Summary (mean +/- SD) per metric x category ──────────────
summary_df <- df %>%
  group_by(metric, category) %>%
  summarise(mean_val = mean(value), sd_val = sd(value), n = n(), .groups = "drop")

cat("\nSummary (mean +/- SD):\n")
print(summary_df)

df_display <- select_representative_points(df)
cat("\nPoints shown on plot per group (representative subsample, minimum 10):\n")
print(table(df_display$metric, df_display$category))

# ── 3. Statistical test per metric: conditional on normality ────
run_conditional_test <- function(data, label) {
  core_vals  <- filter(data, category == "Core")$value
  shell_vals <- filter(data, category == "Shell")$value

  shapiro_core  <- shapiro.test(core_vals)
  shapiro_shell <- shapiro.test(shell_vals)

  cat("\n[", label, "] Shapiro-Wilk: Core p =", signif(shapiro_core$p.value, 3),
      "| Shell p =", signif(shapiro_shell$p.value, 3), "\n")

  both_normal <- shapiro_core$p.value > 0.05 && shapiro_shell$p.value > 0.05

  if (both_normal) {
    test_result <- t.test(value ~ category, data = data)
    test_name <- "Welch's t-test"
  } else {
    test_result <- wilcox.test(value ~ category, data = data, exact = FALSE)
    test_name <- "Wilcoxon rank-sum test"
  }

  p_val <- test_result$p.value
  cat("[", label, "] Using", test_name, "-> p =", signif(p_val, 4),
      " (", format_pval(p_val), ")\n")

  tibble(metric = label, p_val = p_val, p_label = format_pval(p_val))
}

stats_df <- bind_rows(
  run_conditional_test(filter(df, metric == "dN"),      "dN"),
  run_conditional_test(filter(df, metric == "dS"),      "dS"),
  run_conditional_test(filter(df, metric == "dN/dS"),   "dN/dS")
) %>%
  mutate(metric = factor(metric, levels = c("dN", "dS", "dN/dS")))

# ── 4. Per-panel bracket height ──────────────────────────────────
panel_max <- bind_rows(
  summary_df %>% mutate(top = mean_val + sd_val) %>% select(metric, top),
  df_display %>% rename(top = value) %>% select(metric, top)
) %>%
  group_by(metric) %>%
  summarise(max_y = max(top, na.rm = TRUE), .groups = "drop")

stats_df <- stats_df %>%
  left_join(panel_max, by = "metric") %>%
  mutate(
    bracket_y   = max_y * 1.12,
    bracket_gap = max_y * 0.05
  )

# ── 5. Plot: thin dodged bars, black-and-white points, error bars ──
# drawn LAST so they render frontmost (fix applied per plot_heaps_law.R)
dodge_width <- 0.45
bar_width   <- 0.38
cap_width   <- 0.12

summary_df <- summary_df %>% mutate(x_dummy = 1)
df_display <- df_display %>% mutate(x_dummy = 1)

p <- ggplot() +
  # Bar = mean -- flat fill, NO outline
  geom_col(data = summary_df, aes(x = x_dummy, y = mean_val, fill = category),
           position = position_dodge(width = dodge_width),
           width = bar_width, color = NA, alpha = 0.9) +
  # Representative points: black-and-white, shape by category
  geom_point(data = df_display, aes(x = x_dummy, y = value, shape = category, fill = category),
             position = position_jitterdodge(dodge.width = dodge_width,
                                              jitter.width = 0.05, jitter.height = 0),
             color = "black", size = 1.8, stroke = 0.6) +
  # Capped SD error bar ("I" shape), floored at 0, drawn LAST (frontmost)
  geom_errorbar(data = summary_df,
                aes(x = x_dummy, ymin = pmax(mean_val - sd_val, 0), ymax = mean_val + sd_val,
                    group = category),
                position = position_dodge(width = dodge_width),
                width = cap_width, linewidth = 0.6, color = "black") +
  # Significance: actual p-value, plain line spanning the two dodge slots
  geom_segment(data = stats_df,
               aes(x = 1 - dodge_width / 4, xend = 1 + dodge_width / 4,
                   y = bracket_y, yend = bracket_y),
               inherit.aes = FALSE, linewidth = 0.5, color = "black") +
  geom_text(data = stats_df, aes(x = 1, y = bracket_y + bracket_gap, label = p_label),
            inherit.aes = FALSE, size = 3.3, family = "Helvetica") +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = category_palette) +
  scale_shape_manual(values = category_shape) +
  scale_x_continuous(limits = c(1 - dodge_width, 1 + dodge_width)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = NULL, y = "Evolutionary Rate") +
  theme_figure() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = c(1, 0), color = "black")))

# ── 6. Save ────────────────────────────────────────────────────
fig_width_in  <- 7
fig_height_in <- 5.5

out_path <- file.path(base_dir, "dnds_core_shell.png")
ggsave(out_path, p, width = fig_width_in, height = fig_height_in, dpi = 300, bg = "white")

cat("\nSaved:", out_path, "\n")
p
