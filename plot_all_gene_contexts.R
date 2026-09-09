# ============================================================
# All 6 LyasA/LyasB genome-specific gene contexts, combined into ONE
# figure, stacked as tracks sharing a COMMON coordinate scale (same
# bp-per-inch across every panel) -- standard genome-browser
# convention. This fixes the earlier "stretched" artifact: forcing
# every panel to the same fixed WIDTH regardless of actual bp span
# made tiny regions (791bp) look as visually "full" as huge ones
# (20,000bp). Now every panel shares the same x-axis SPAN; small
# regions simply show as a short filled stretch on an otherwise-blank
# shared axis, rather than being artificially stretched to fill it.
#
# Legend kept compact (short gene IDs only); full product descriptions
# printed + saved separately per region, not crammed into the figure.
# ============================================================

library(tidyverse)
library(colorspace)
library(jsonlite)
library(patchwork)

context_data <- fromJSON(path.expand("~/Downloads/all_gene_contexts.json"), simplifyVector = FALSE)

make_arrow <- function(gene_id, start, end, strand, half_h, head_frac = 0.22) {
  len <- end - start
  head_len <- max(len * head_frac, 1)
  body_h <- half_h * 0.55
  if (strand == "+") {
    tibble(x = c(start, end - head_len, end - head_len, end, end - head_len, end - head_len, start),
           y = c(body_h, body_h, half_h, 0, -half_h, -body_h, -body_h), gene_id = gene_id)
  } else {
    tibble(x = c(end, start + head_len, start + head_len, start, start + head_len, start + head_len, end),
           y = c(body_h, body_h, half_h, 0, -half_h, -body_h, -body_h), gene_id = gene_id)
  }
}

# ── Shared scale: every panel spans the SAME width (largest region's
# span), so bp-per-inch is identical across all 6 -- the actual fix
# for the stretching problem. ──────────────────────────────────────
all_spans <- sapply(context_data, function(e) e$window_end - e$window_start)
MAX_SPAN <- max(all_spans)
cat("Region spans (bp):\n"); print(all_spans)
cat("Shared axis span used for all panels:", MAX_SPAN, "bp\n\n")

theme_track <- theme_minimal(base_family = "Helvetica") +
  theme(
    panel.grid    = element_blank(),
    axis.line.y   = element_blank(), axis.text.y = element_blank(),
    axis.ticks.y  = element_blank(), axis.title.y = element_blank(),
    axis.line.x   = element_line(color = "black", linewidth = 0.3),
    axis.ticks.x  = element_line(color = "black", linewidth = 0.25),
    axis.text.x   = element_text(size = 7.5, color = "black"),
    axis.title.x  = element_blank(),
    plot.title    = element_text(face = "bold", size = 9, hjust = 0),
    legend.position = "right",
    legend.title      = element_text(face = "bold", size = 7.5),
    legend.text        = element_text(size = 6.2),
    legend.key.size      = unit(6, "pt"),
    plot.margin   = margin(2, 10, 2, 10)
  )

build_track <- function(entry, region_name) {
  target_ids <- unlist(entry$target_genes)

  genes <- map_dfr(entry$genes, function(g) {
    tibble(gene_id = g$gene_id, start = g$start, end = g$end,
           strand = g$strand, product = g$product)
  }) %>%
    mutate(is_target = gene_id %in% target_ids,
           legend_label = paste0(ifelse(is_target, "\u2605 ", ""), gene_id))

  n_genes <- nrow(genes)
  gene_colors <- setNames(qualitative_hcl(n_genes, palette = "Dark 3"), genes$legend_label)

  window_start <- entry$window_start
  window_end <- window_start + MAX_SPAN   # SAME span for every panel
  half_h <- MAX_SPAN * 0.008 / 18

  arrows <- pmap_dfr(genes, function(gene_id, start, end, strand, legend_label, ...)
    make_arrow(gene_id, start, end, strand, half_h)) %>%
    left_join(genes %>% select(gene_id, legend_label), by = "gene_id")

  legend_ncol <- max(1, ceiling(n_genes / 10))

  ggplot() +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.25) +
    geom_polygon(data = arrows, aes(x = x, y = y, group = gene_id, fill = legend_label),
                 color = "white", linewidth = 0.3) +
    scale_fill_manual(values = gene_colors, name = NULL) +
    scale_x_continuous(limits = c(entry$window_start, window_end),
                        breaks = scales::breaks_pretty(n = 4), labels = scales::comma,
                        expand = c(0, 0)) +
    coord_cartesian(ylim = c(-half_h * 2.5, half_h * 2.5), clip = "off") +
    guides(fill = guide_legend(ncol = legend_ncol, byrow = TRUE)) +
    labs(title = paste0(gsub("_", " ", region_name), "  (", entry$contig, ", ",
                          format(entry$contig_length, big.mark=","), "bp total)")) +
    theme_track
}

tracks <- imap(context_data, build_track)

combined <- wrap_plots(tracks, ncol = 1) +
  plot_annotation(
    title = "LyasA / LyasB genome-specific gene regions -- shared coordinate scale",
    subtitle = paste0("All panels span ", format(MAX_SPAN, big.mark = ","),
                        "bp (largest region) -- smaller regions show as a shorter filled stretch, not stretched to fill the panel"),
    theme = theme(plot.title = element_text(face = "bold", size = 13, family = "Helvetica"),
                   plot.subtitle = element_text(size = 9.5, color = "grey40", family = "Helvetica"))
  )

out_path <- path.expand("~/Downloads/all_gene_contexts_combined.png")
ggsave(out_path, combined, width = 12, height = (2.2 / 3) * length(tracks) + 1, dpi = 300, bg = "white")
cat("Saved:", out_path, "\n")

# Full reference tables, same as before -- separate from the compact
# on-figure legends
for (region in names(context_data)) {
  entry <- context_data[[region]]
  target_ids <- unlist(entry$target_genes)
  ref_table <- map_dfr(entry$genes, function(g) {
    tibble(gene_id = g$gene_id, is_target = g$gene_id %in% target_ids, product = g$product)
  })
  write_csv(ref_table, path.expand(paste0("~/Downloads/", region, "_gene_reference.csv")))
}
cat("Gene reference CSVs saved per region.\n")
