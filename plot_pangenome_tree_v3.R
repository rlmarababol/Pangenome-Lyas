# ============================================================
# Pangenome phylogenomic tree + gene presence/absence heatmap
# L. yasudae core genome tree (FastTree) + Roary-style prevalence bins
#
# v3 changes: higher-resolution raster output + a vector PDF export
# alongside it, per request. Per the style guide (Section 9), vector
# output is preferred for lines/points/text, with raster reserved for
# panels with very large tile counts -- this matrix has ~5,300+ tile
# columns, so the PDF here may render slowly to open/print and can be
# a large file; the PNG remains the more practical file for quick
# viewing/sharing, the PDF for print-quality/editable-text use.
#
# Revision notes carried over from v2:
#   - Genome tip labels right-aligned flush against the heatmap edge
#     (geom_tiplab(align = TRUE)), with extra horizontal room so long
#     names are never clipped.
#   - Matrix tiles show PRESENCE/ABSENCE only (no partition color-coding
#     in the tile fill): Absent = faded violet, Present = dark violet --
#     a tint ladder of the guide's Lavender/Purple anchor (Section 2.2),
#     since present/absent is one entity at two intensities, not two
#     unrelated categories.
#   - The four conventional Roary prevalence tiers -- Core, Soft core,
#     Shell, Cloud -- are shown as facet labels above column blocks
#     (Section 8: "facet by the natural grouping variable"), not as
#     tile colors. Thresholds: Core >=99%, Soft core 95-99%,
#     Shell 15-95%, Cloud <15% of genomes.
#   - Cell shape: with ~5,300+ family columns against only 10 genome
#     rows, literal 1:1 square cells would need an image many feet
#     wide (not practically viewable). Cells are pulled meaningfully
#     closer to square via a wider matrix panel PLUS thickened white
#     horizontal separators between genome rows (visually "eating into"
#     each row's apparent height), rather than resizing alone.
# ============================================================

library(tidyverse)
library(ggtree)
library(treeio)
library(patchwork)
library(colorspace)

# ── 0. Shared palette + theme (single source of truth) ────────
# Presence/absence is really one entity at two "intensities" (gene family
# detected vs. not) -- so per the guide's Section 2.2 (ordinal/dose-like
# variables get a TINT LADDER of one hue, not two unrelated hues), both
# colors are built from the guide's Lavender/Purple anchor (#B491E5):
# full-strength dark violet for Present, a faded light tint for Absent.
base_violet  <- "#B491E5"                                        # guide anchor
present_col  <- colorspace::darken(base_violet,  amount = 0.45)  # dark violet
absent_col   <- colorspace::lighten(base_violet, amount = 0.55)  # faded violet

presence_palette <- c(Absent = absent_col, Present = present_col)

theme_figure <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "Helvetica") +
    theme(
      panel.grid        = element_blank(),
      axis.line         = element_line(color = "black", linewidth = 0.4),
      axis.ticks        = element_line(color = "black", linewidth = 0.3),
      axis.ticks.length = unit(2, "pt"),
      axis.title        = element_text(face = "bold", size = base_size),
      axis.text         = element_text(size = base_size - 2, color = "black"),
      legend.title       = element_text(face = "bold", size = base_size - 1),
      legend.text        = element_text(size = base_size - 2),
      legend.background  = element_blank(),
      legend.key         = element_blank(),
      strip.background   = element_blank(),
      strip.text         = element_text(face = "bold", size = base_size - 1),
      plot.tag            = element_text(face = "bold", size = base_size + 4),
      plot.margin         = margin(6, 18, 6, 6)
    )
}

# ── 1. Load data ─────────────────────────────────────────────
base_dir  <- "~/Downloads/macrogen_1c/pangenome/pangenome_results"
tree_file <- file.path(base_dir, "core_genome.tree")
rtab_file <- file.path(base_dir, "matrix_export/gene_presence_absence.Rtab")

tree <- read.tree(tree_file)
rtab <- read_tsv(rtab_file, show_col_types = FALSE)
gene_col <- colnames(rtab)[1]  # usually "Gene" or similar family-ID column

genome_cols <- setdiff(colnames(rtab), gene_col)
n_genomes   <- length(genome_cols)

cat("Tree tips:", length(tree$tip.label), "\n")
cat("Rtab:", nrow(rtab), "families x", n_genomes, "genomes\n")

# ── 2. Classify each family into standard Roary prevalence tiers ──
# Core: >=99% | Soft core: 95-99% | Shell: 15-95% | Cloud: <15%
# NOTE: with only 10 genomes, prevalence only takes discrete 10% steps,
# so the 95-99% "Soft core" band may come out empty/very sparse -- that
# is a property of the discrete genome count, not a script bug. It will
# fill in naturally as more genomes are added to the pangenome.
category_levels <- c("Core", "Soft core", "Shell", "Cloud")

mat <- rtab %>%
  rename(family = !!gene_col) %>%
  mutate(
    prevalence     = rowSums(across(all_of(genome_cols))),
    prevalence_pct = 100 * prevalence / n_genomes,
    category = case_when(
      prevalence_pct >= 99 ~ "Core",
      prevalence_pct >= 95 ~ "Soft core",
      prevalence_pct >= 15 ~ "Shell",
      TRUE                 ~ "Cloud"
    ),
    category = factor(category, levels = category_levels)
  ) %>%
  arrange(category, desc(prevalence))

cat("\nFamilies per tier:\n")
print(table(mat$category))

mat_long <- mat %>%
  select(family, category, all_of(genome_cols)) %>%
  pivot_longer(all_of(genome_cols), names_to = "genome", values_to = "presence") %>%
  mutate(
    family   = factor(family, levels = unique(mat$family)),
    presence = factor(presence, levels = c(0, 1), labels = c("Absent", "Present"))
  )

# ── 3. Tree panel (tag "a") ───────────────────────────────────
# align = TRUE right-aligns every tip label at a common x position
# (dotted connector for tips shorter than the deepest branch) so names
# read as one clean flush column right next to the heatmap. hexpand()
# reserves extra right-margin room so long names never get clipped.
# If your ggtree version lacks hexpand(), replace that line with:
#   + xlim(0, max(p_tree$data$x, na.rm = TRUE) * 1.4)
p_tree <- ggtree(tree, linewidth = 0.5, color = "black") +
  geom_tiplab(align = TRUE, linetype = "dotted", linewidth = 0.25,
              size = 3.2, family = "Helvetica", offset = 0.001) +
  hexpand(0.35, direction = 1) +
  theme_figure() +
  theme(
    axis.line.y  = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(x = "Substitutions / site")

# Lock genome order so heatmap rows stay pixel-aligned with tree tips.
tip_order <- get_taxa_name(p_tree)
mat_long <- mat_long %>% mutate(genome = factor(genome, levels = rev(tip_order)))

# ── 4. Heatmap panel (tag "b") ────────────────────────────────
# Plain presence/absence fill only. Core/Soft core/Shell/Cloud appear
# as facet strip labels, with each block's width proportional to its
# family count (space = "free_x") rather than encoded via tile color.
#
# CELL-SHAPE FIX: geom_tile()'s own `color`/`linewidth` draws a border
# on all four edges at once, so it can't selectively thicken just the
# row separators. Instead, tile borders are turned off (color = NA) and
# thick white horizontal lines are drawn explicitly at each row
# boundary -- this "eats into" the visible height of every row without
# touching column width, which is what actually pulls the apparent
# cell shape toward square when there are far more columns than rows.
row_boundaries <- tibble(y = seq(1.5, n_genomes - 0.5, by = 1))

p_matrix <- ggplot(mat_long, aes(x = family, y = genome, fill = presence)) +
  geom_tile(color = NA) +
  geom_hline(data = row_boundaries, aes(yintercept = y),
             color = "white", linewidth = 1.3, inherit.aes = FALSE) +
  facet_grid(~ category, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = presence_palette, name = "Gene family") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_figure() +
  theme(
    axis.line       = element_blank(),
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    axis.title      = element_blank(),
    panel.spacing   = unit(4, "pt"),
    legend.position = "right"
  )

# ── 5. Compose: aligned panels, bold lowercase tags ───────────
combined <- (p_tree | p_matrix) +
  plot_layout(widths = c(1, 4)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 15, family = "Helvetica"))

# ── 6. Save: high-resolution PNG + vector PDF ─────────────────
# fig_width_in / fig_height_in control both outputs identically, so PNG
# and PDF always show the exact same layout/proportions.
fig_width_in  <- 40
fig_height_in <- 6

# PNG: raster, higher DPI than before for print-quality zoom/crop.
# 600 dpi is roughly print-shop quality; drop to 300-400 if the file
# gets unwieldy to open/share.
png_path <- file.path(base_dir, "phylogenomic_tree_pangenome.png")
ggsave(png_path, combined, width = fig_width_in, height = fig_height_in,
       dpi = 600, bg = "white")
cat("\nSaved PNG:", png_path, "\n")

# PDF: vector, via cairo_pdf for correct font embedding (Section 9:
# vector output for lines/points/text). With ~5,300+ tile columns this
# file may be large and slow to open in Illustrator/Preview -- that is
# the expected tradeoff for a heatmap this dense rendered as true
# vector rather than a raster exception.
pdf_path <- file.path(base_dir, "phylogenomic_tree_pangenome.pdf")
ggsave(pdf_path, combined, width = fig_width_in, height = fig_height_in,
       device = cairo_pdf, bg = "white")
cat("Saved PDF:", pdf_path, "\n")

combined
