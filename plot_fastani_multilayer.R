# ============================================================
# Multi-layer whole-genome synteny, all 10 genomes, tree order
# (LyasH top, LyasI bottom) -- using fastANI's ACTUAL reciprocal
# mapping data (--visualize / .visual files, confirmed correct),
# rendered via gggenomes' geom_link() for smooth curved ribbons
# (genoPlotR, used for the single-pair version, draws straight-edged
# polygons only -- no native curve/smoothing option).
#
# KEY: fastANI's .visual coordinates are CONCATENATED across all
# contigs within each genome (confirmed via a FastANI GitHub issue
# showing this exact behavior) -- so each genome becomes ONE clean
# backbone "pseudo-sequence" here, sidestepping the fragmented-contig
# mess (25-89 contigs per genome) that caused problems in every
# earlier attempt.
#
# .visual format (BLAST-tabular-compatible, confirmed via FastANI
# issue #133): query, subject, pident, NA, NA, NA, qstart, qend,
# sstart, send, NA, NA
# ============================================================

library(gggenomes)
library(tidyverse)

syn_dir <- "~/Downloads/macrogen_1c/pangenome/pangenome_results/synteny"
visual_dir <- file.path(syn_dir, "fastani_visual")
genome_dir <- file.path(syn_dir, "ani_genomes")

TREE_ORDER <- c("LyasH", "LyasD", "LyasG", "LyasE", "LyasC",
                 "LyasJ", "LyasF", "LyasB", "LyasA", "LyasI")
PAIRS <- list(c("LyasH","LyasD"), c("LyasD","LyasG"), c("LyasG","LyasE"),
              c("LyasE","LyasC"), c("LyasC","LyasJ"), c("LyasJ","LyasF"),
              c("LyasF","LyasB"), c("LyasB","LyasA"), c("LyasA","LyasI"))

# ── 1. Total (concatenated) length per genome, from .fai files ───
genome_lengths <- map_dfr(TREE_ORDER, function(g) {
  fai <- read_tsv(file.path(genome_dir, paste0(g, ".fasta.fai")),
                   col_names = FALSE, show_col_types = FALSE)
  tibble(bin_id = g, length = sum(fai$X2))
})

cat("Genome (concatenated) lengths:\n")
print(genome_lengths)

seqs <- genome_lengths %>%
  mutate(seq_id = bin_id, bin_id = factor(bin_id, levels = TREE_ORDER)) %>%
  arrange(bin_id) %>%
  select(seq_id, bin_id, length)

# ── 2. Parse all 9 .visual files (query=g1, reference=g2 per pair) ──
parse_visual <- function(g1, g2) {
  path <- file.path(visual_dir, paste0(g1, "_vs_", g2, ".out.visual"))
  df <- read_tsv(path, col_names = c("query", "subject", "pident", "na1", "na2", "na3",
                                       "qstart", "qend", "sstart", "send", "na4", "na5"),
                  show_col_types = FALSE)
  df %>%
    transmute(seq_id = g1, start = qstart, end = qend,
              seq_id2 = g2, start2 = sstart, end2 = send,
              pident = as.numeric(pident))
}

all_links <- map2_dfr(map_chr(PAIRS, 1), map_chr(PAIRS, 2), parse_visual)

cat("\nTotal mapping segments across all 9 pairs:", nrow(all_links), "\n")
cat("%identity summary:\n")
print(summary(all_links$pident))

# ── 3. Plot with smooth curved ribbons ────────────────────────────
gg <- gggenomes(seqs = seqs, links = all_links)

p <- gg +
  geom_seq(linewidth = 1, color = "black") +
  geom_link_curved(aes(fill = pident), color = NA, alpha = 0.75) +
  scale_fill_gradient(low = "#EAF1FB", high = "#66A8E5", name = "% identity") +
  theme_minimal(base_size = 11, base_family = "Helvetica") +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    legend.position = "right"
  )

out_path <- file.path(syn_dir, "fastani_multilayer_synteny.png")
ggsave(out_path, p, width = 14, height = 8, dpi = 300, bg = "white")

cat("\nSaved:", out_path, "\n")
p
