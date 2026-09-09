# ============================================================
# Pairwise ANI heatmap, 10 L. yasudae genomes
# Row/column order matches plot_synteny.R exactly (H at top, I at
# bottom) so the two panels align if placed together.
#
# fastANI's raw output is NOT perfectly symmetric (A-vs-B can differ
# slightly from B-vs-A, since it fragments the QUERY genome, not the
# reference) -- both directions are averaged per pair for a clean
# symmetric matrix, the standard approach for this kind of heatmap.
# ============================================================

library(tidyverse)

base_dir <- "~/Downloads/macrogen_1c/pangenome/pangenome_results/synteny"

TREE_ORDER <- c("LyasH", "LyasD", "LyasG", "LyasE", "LyasC",
                 "LyasJ", "LyasF", "LyasB", "LyasA", "LyasI")

theme_figure <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "Helvetica") +
    theme(
      panel.grid  = element_blank(),
      axis.title  = element_blank(),
      axis.text   = element_text(size = base_size - 1, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(face = "bold", size = base_size - 1),
      legend.text  = element_text(size = base_size - 2)
    )
}

# ── 1. Load raw fastANI output, extract clean genome codes ──────
raw <- read_tsv(file.path(base_dir, "ani_all_vs_all.tsv"),
                 col_names = c("query", "ref", "ani", "matched", "total"),
                 show_col_types = FALSE)

extract_code <- function(path) {
  str_remove(basename(path), "\\.fasta$")
}

df <- raw %>%
  mutate(genome1 = extract_code(query), genome2 = extract_code(ref)) %>%
  select(genome1, genome2, ani)

cat("Rows loaded:", nrow(df), "\n")
cat("Unique genomes seen:", n_distinct(c(df$genome1, df$genome2)), "\n")

# ── 2. Symmetrize: average both directions per unordered pair ───
symmetric <- df %>%
  mutate(
    g_lo = pmin(genome1, genome2),
    g_hi = pmax(genome1, genome2)
  ) %>%
  group_by(g_lo, g_hi) %>%
  summarise(ani_mean = mean(ani), n_dir = n(), .groups = "drop")

cat("\nUnordered pairs found:", nrow(symmetric),
    "(expect 10 self + 45 off-diagonal = 55)\n")
cat("Pairs with only 1 direction (check if any pair failed fastANI's threshold):\n")
print(symmetric %>% filter(n_dir < 2, g_lo != g_hi))

# ── 3. Expand to full matrix (both triangles + diagonal) ────────
full_matrix <- bind_rows(
  symmetric %>% rename(genome1 = g_lo, genome2 = g_hi),
  symmetric %>% rename(genome1 = g_hi, genome2 = g_lo) %>% filter(genome1 != genome2)
) %>%
  mutate(
    genome1 = factor(genome1, levels = rev(TREE_ORDER)),  # rev() -> H renders at top
    genome2 = factor(genome2, levels = TREE_ORDER)         # left-to-right, H first
  )

# ── 4. Plot ────────────────────────────────────────────────────
p <- ggplot(full_matrix, aes(x = genome2, y = genome1, fill = ani_mean)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", ani_mean)), size = 2.6, color = "black") +
  scale_fill_gradient(low = "#F5F0FC", high = "#B491E5", name = "ANI (%)") +
  scale_x_discrete(position = "top") +
  coord_fixed() +
  theme_figure()

out_path <- file.path(base_dir, "ani_heatmap.png")
ggsave(out_path, p, width = 7, height = 6.5, dpi = 300, bg = "white")

cat("\nSaved:", out_path, "\n")
p
