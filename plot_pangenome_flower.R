# ============================================================
# Pangenome flower plot: center = Core (shared across all 10 genomes),
# petals = each of the 9 non-reference genomes' dispensable/unique
# gene count. NOT a true Venn diagram (no pairwise intersection
# regions) -- matches the simpler standard pangenomics flower
# convention the user referenced, with real elongated teardrop petal
# shapes (ggforce::geom_ellipse, rotated/offset per genome) rather
# than simple point markers.
#
# 9 petals (LyasA-H, J) -- LyasI excluded as the fixed reference
# genome throughout this project, consistent with the PCA color
# legend this palette is recovered from.
# ============================================================

library(tidyverse)
library(ggforce)
library(colorspace)

base_dir <- "~/Downloads/macrogen_1c/pangenome/pangenome_results"

genome_codes_9 <- sort(c("LyasA","LyasB","LyasC","LyasD","LyasE","LyasF","LyasG","LyasH","LyasJ"))
genome_palette <- setNames(qualitative_hcl(9, palette = "Set 2"), genome_codes_9)

# ── Load data ──────────────────────────────────────────────────────
rtab <- read_tsv(file.path(base_dir, "matrix_export/gene_presence_absence.Rtab"), show_col_types = FALSE)
gene_col <- colnames(rtab)[1]
genome_cols <- setdiff(colnames(rtab), gene_col)
exact_core_ids <- readLines(file.path(base_dir, "matrix_export/partitions/exact_core.txt"))
exact_core_count <- length(exact_core_ids)

fam_data <- rtab %>%
  rename(family = !!gene_col) %>%
  mutate(prevalence = rowSums(across(all_of(genome_cols))))

unique_counts <- fam_data %>%
  filter(prevalence == 1) %>%
  rowwise() %>%
  mutate(genome = genome_cols[which(c_across(all_of(genome_cols)) == 1)]) %>%
  ungroup() %>%
  count(genome, name = "unique_genes") %>%
  filter(genome %in% genome_codes_9) %>%   # 9 petals only, LyasI excluded
  right_join(tibble(genome = genome_codes_9), by = "genome") %>%
  mutate(unique_genes = replace_na(unique_genes, 0)) %>%
  arrange(match(genome, genome_codes_9))

cat("Core (center):", exact_core_count, "\n")
cat("Genome-specific (petal) counts:\n")
print(unique_counts)

# ── Petal geometry: ellipse per genome, offset + rotated outward ──
n <- nrow(unique_counts)
core_radius <- 2.2
petal_offset <- 2.6      # distance from center to each petal's center
petal_length <- 2.4      # semi-major axis (elongation, outward)
petal_width  <- 1.15     # semi-minor axis (petal thickness)

petals <- unique_counts %>%
  mutate(
    angle_deg = seq(90, 90 - 360 * (n - 1) / n, length.out = n),
    angle_rad = angle_deg * pi / 180,
    x0 = petal_offset * cos(angle_rad),
    y0 = petal_offset * sin(angle_rad),
    label_x = (petal_offset + petal_length + 0.5) * cos(angle_rad),
    label_y = (petal_offset + petal_length + 0.5) * sin(angle_rad),
    count_label_x = (petal_offset + 0.5) * cos(angle_rad),
    count_label_y = (petal_offset + 0.5) * sin(angle_rad)
  )

theme_flower <- theme_void(base_family = "Helvetica") +
  theme(legend.position = "none", plot.margin = margin(25, 25, 25, 25))

p <- ggplot() +
  # petals (drawn first, so the central circle sits on top of their bases)
  geom_ellipse(data = petals,
               aes(x0 = x0, y0 = y0, a = petal_length, b = petal_width,
                   angle = angle_rad, fill = genome),
               color = "white", linewidth = 0.6, alpha = 0.85) +
  scale_fill_manual(values = genome_palette) +
  # petal count labels
  geom_text(data = petals, aes(x = count_label_x, y = count_label_y, label = unique_genes),
            size = 3.6, family = "Helvetica", fontface = "bold", color = "black") +
  # genome name labels, outside the petals
  geom_text(data = petals, aes(x = label_x, y = label_y, label = genome),
            size = 3.8, family = "Helvetica", fontface = "bold") +
  # central core circle, drawn on top
  ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = core_radius), fill = "#B491E5", color = "white", linewidth = 0.8) +
  geom_text(aes(x = 0, y = 0, label = paste0("Core\n", format(exact_core_count, big.mark = ","))),
            size = 5, family = "Helvetica", fontface = "bold", color = "white", lineheight = 0.9) +
  coord_fixed(clip = "off", xlim = c(-7, 7), ylim = c(-7, 7)) +
  theme_flower

out_path <- file.path(base_dir, "pangenome_flower.png")
ggsave(out_path, p, width = 10, height = 10, dpi = 300, bg = "white")
cat("\nSaved:", out_path, "\n")
