# ============================================================
# Top 10 gene families by total gene count, per partition category --
# now as THREE SEPARATE standalone figures (Core, Shell, Cloud) instead
# of one faceted combined plot.
#
# Same data/logic as the original combined version: "number of genes"
# = total individual gene copies belonging to that family across ALL
# 10 genomes, labeled with product annotations (not raw locus tags).
# ============================================================

library(tidyverse)
library(tidytext)

base_dir <- "~/Downloads/macrogen_1c/pangenome/pangenome_results"
syn_dir  <- file.path(base_dir, "synteny")

theme_figure <- function(base_size = 12) {
  theme_minimal(base_size = base_size, base_family = "Helvetica") +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor    = element_blank(),
      axis.line.x        = element_line(color = "black", linewidth = 0.4),
      axis.ticks.x       = element_line(color = "black", linewidth = 0.3),
      axis.line.y         = element_blank(),
      axis.ticks.y        = element_blank(),
      axis.title           = element_text(face = "bold", size = base_size),
      axis.text             = element_text(size = base_size - 1.5, color = "black"),
      plot.title             = element_text(face = "bold", size = base_size + 2),
      legend.position          = "none"
    )
}

category_palette <- c("Core" = "#B491E5", "Shell" = "#66A8E5", "Cloud" = "#E4749B")

# ── Load genes, product annotations, partition family lists ──────
genes <- read_csv(file.path(syn_dir, "genes.tsv"), show_col_types = FALSE)
products <- read_csv(file.path(syn_dir, "gene_products.csv"), show_col_types = FALSE)

partitions_dir <- file.path(base_dir, "matrix_export/partitions")
core_ids  <- readLines(file.path(partitions_dir, "persistent.txt"))
shell_ids <- readLines(file.path(partitions_dir, "shell.txt"))
cloud_ids <- readLines(file.path(partitions_dir, "cloud.txt"))

category_map <- bind_rows(
  tibble(family_id = core_ids,  category = "Core"),
  tibble(family_id = shell_ids, category = "Shell"),
  tibble(family_id = cloud_ids, category = "Cloud")
)

gene_counts <- genes %>%
  filter(!is.na(family_id)) %>%
  count(family_id, name = "n_genes") %>%
  inner_join(category_map, by = "family_id")

family_products <- genes %>%
  filter(!is.na(family_id)) %>%
  left_join(products, by = "gene_id") %>%
  filter(!is.na(product)) %>%
  count(family_id, product, sort = TRUE) %>%
  group_by(family_id) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(family_id, product)

gene_counts <- gene_counts %>%
  left_join(family_products, by = "family_id") %>%
  mutate(
    product = replace_na(product, "hypothetical protein"),
    label = paste0(product, " (", family_id, ")")
  )

top10 <- gene_counts %>%
  group_by(category) %>%
  slice_max(n_genes, n = 10, with_ties = FALSE) %>%
  ungroup()

# ── Render one standalone plot per category ───────────────────────
for (cat in c("Core", "Shell", "Cloud")) {
  df <- top10 %>%
    filter(category == cat) %>%
    mutate(label = fct_reorder(label, n_genes))

  p <- ggplot(df, aes(x = n_genes, y = label)) +
    geom_col(width = 0.7, color = NA, fill = category_palette[cat], alpha = 0.9) +
    labs(x = "Number of genes", y = NULL, title = paste0("Top 10 ", cat, " gene families")) +
    theme_figure()

  out_path <- file.path(base_dir, paste0("top10_", tolower(cat), "_gene_families.png"))
  ggsave(out_path, p, width = 9, height = 5, dpi = 300, bg = "white")
  cat_msg <- paste0("Saved: ", out_path, "\n")
  message(cat_msg)
}

cat("\nAll three saved separately.\n")
