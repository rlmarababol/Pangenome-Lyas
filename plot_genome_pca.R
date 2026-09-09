# ============================================================
# PCA of each genome's evolutionary profile, separately for its
# Core-gene metrics (circle) and Shell-gene metrics (square)
#
# Each genome contributes TWO observations to the PCA: one vector of
# {SNP/kb, Indel/kb, dN, dS, dN/dS, Frameshift/kb, In-frame-indel/kb}
# computed from its Core genes, and a second vector of the same seven
# metrics computed from its Shell genes (both already validated in
# earlier scripts -- no new computation, just reshaping existing
# per-genome results into one feature matrix).
#
# Genome identity needs its own hue (9 genomes) -- more categories than
# the style guide's small anchor palette (~6) is designed for, so a
# proper qualitative palette is used here deliberately, not as a
# drift from the system. Category (Core/Shell) is the sub-grouping,
# encoded by shape per Section 1.4's convention.
#
# RELIABILITY FILTER: same convention used throughout this pipeline --
# genome x category combinations backed by too few compared families
# give unreliable metrics and are excluded (affects LyasD/G Shell,
# LyasH Shell entirely absent -- the fragmented-assembly trio).
# ============================================================

library(tidyverse)

MIN_FAMILIES <- 10
base_dir <- "~/Downloads/macrogen_1c/pangenome/pangenome_results"

theme_figure <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "Helvetica") +
    theme(
      panel.grid        = element_blank(),
      axis.line         = element_line(color = "black", linewidth = 0.4),
      axis.ticks        = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(2, "pt"),
      axis.title        = element_text(face = "bold", size = base_size),
      axis.text         = element_text(size = base_size - 1, color = "black"),
      legend.title       = element_text(face = "bold", size = base_size - 1),
      legend.text        = element_text(size = base_size - 2),
      legend.background  = element_blank(),
      legend.key          = element_blank()
    )
}

# ── 1. Load and merge per-genome metrics already computed ───────
snp_indel <- read_csv(file.path(base_dir, "snp_indel_per_genome.csv"), show_col_types = FALSE) %>%
  select(genome, category, n_families_compared, total_aln_len_bp, snp_per_kb, indel_per_kb)

dnds <- read_csv(file.path(base_dir, "dnds_per_genome.csv"), show_col_types = FALSE) %>%
  select(genome, category, dN, dS, dN_dS)

large_effect <- read_csv(file.path(base_dir, "large_effect_summary.csv"), show_col_types = FALSE) %>%
  select(genome, category, n_frameshift, n_inframe_indel)

merged <- snp_indel %>%
  left_join(dnds, by = c("genome", "category")) %>%
  left_join(large_effect, by = c("genome", "category")) %>%
  mutate(
    category = factor(category, levels = c("Core", "Shell")),
    frameshift_per_kb    = n_frameshift / (total_aln_len_bp / 1000),
    inframe_indel_per_kb = n_inframe_indel / (total_aln_len_bp / 1000)
  )

cat("All genome x category rows:\n")
print(merged %>% select(genome, category, n_families_compared, snp_per_kb, indel_per_kb,
                          dN, dS, dN_dS, frameshift_per_kb, inframe_indel_per_kb))

excluded <- merged %>% filter(n_families_compared < MIN_FAMILIES)
cat("\nExcluded for unreliable low family count (< ", MIN_FAMILIES, "):\n", sep = "")
print(excluded %>% select(genome, category, n_families_compared))

merged <- merged %>% filter(n_families_compared >= MIN_FAMILIES) %>%
  drop_na(snp_per_kb, indel_per_kb, dN, dS, dN_dS, frameshift_per_kb, inframe_indel_per_kb)

cat("\nObservations retained for PCA:", nrow(merged), "\n")

# ── 2. Build feature matrix, run PCA (scaled -- features are on very
# different numeric scales, e.g. dN/dS ~0.1 vs SNP/kb ~10-40) ───────
feature_cols <- c("snp_per_kb", "indel_per_kb", "dN", "dS", "dN_dS",
                   "frameshift_per_kb", "inframe_indel_per_kb")

feature_matrix <- merged %>% select(all_of(feature_cols)) %>% as.matrix()
rownames(feature_matrix) <- paste(merged$genome, merged$category, sep = "_")

pca <- prcomp(feature_matrix, center = TRUE, scale. = TRUE)

var_explained <- summary(pca)$importance[2, ] * 100
cat("\nVariance explained: PC1 =", round(var_explained[1], 1),
    "% | PC2 =", round(var_explained[2], 1), "%\n")

pca_df <- merged %>%
  mutate(PC1 = pca$x[, 1], PC2 = pca$x[, 2])

# ── 3. Plot ───────────────────────────────────────────────────
# Hollow/open markers (color = outline only, no fill) so heavily
# overlapping points stay distinguishable -- matches the reference
# image's fix for the same overlap problem.
library(colorspace)
genome_codes <- sort(unique(pca_df$genome))
genome_palette <- setNames(
  colorspace::qualitative_hcl(length(genome_codes), palette = "Set 2"),
  genome_codes
)
category_shape <- c("Core" = 1, "Shell" = 0)  # open circle, open square

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = genome, shape = category)) +
  # Zero-reference lines, muted and behind the data (Section 6.2)
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  # Group ellipses by category (the "measurement" a genome's profile was
  # built from) -- neutral color so they don't compete with the genome
  # hues; linetype matches the same Category variable already shown by
  # point shape, so the two legends merge into one.
  stat_ellipse(aes(group = category, linetype = category), color = "grey40",
               linewidth = 0.5, level = 0.95, inherit.aes = TRUE) +
  geom_point(size = 3.5, stroke = 1, alpha = 0.9) +
  scale_color_manual(values = genome_palette, name = "Genome") +
  scale_shape_manual(values = category_shape, name = "Category") +
  scale_linetype_manual(values = c("Core" = "solid", "Shell" = "dashed"), name = "Category") +
  labs(
    x = paste0("PC1: ", round(var_explained[1], 0), "% variance"),
    y = paste0("PC2: ", round(var_explained[2], 0), "% variance")
  ) +
  theme_figure()

out_path <- file.path(base_dir, "genome_pca.png")
ggsave(out_path, p, width = 7.5, height = 6, dpi = 300, bg = "white")

cat("\nSaved:", out_path, "\n")
p
