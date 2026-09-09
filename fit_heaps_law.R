# ============================================================
# Heap's Law fit: pangenome_size = kappa * N^gamma
#
# gamma (often called "B" or "alpha" depending on the paper's
# convention -- here explicitly labeled gamma to avoid ambiguity) is
# the parameter that determines openness:
#   gamma close to 1   -> open pangenome, still discovering genes
#                          near-linearly with each new genome
#   gamma close to 0   -> closed pangenome, approaching saturation
#   (by convention, gamma > 0 = open, gamma <= 0 = closed, following
#   Tettelin et al. 2008's original formulation)
#
# Fit uses ALL individual replicate rows (not per-N averages), which
# properly captures sampling variability across the 30 random genome
# orderings, rather than smoothing it away.
# ============================================================

library(tidyverse)

base_dir <- "~/Downloads/macrogen_1c/pangenome/pangenome_results/rarefaction"
rarefaction <- read_csv(file.path(base_dir, "rarefaction.csv"), show_col_types = FALSE)

cat("Total rows (replicate observations):", nrow(rarefaction), "\n")
cat("Genome-count range:", min(rarefaction$genomes_count), "-", max(rarefaction$genomes_count), "\n")
cat("Rows per genome_count:\n")
print(rarefaction %>% count(genomes_count))

# ── Fit Heap's Law: pangenome = kappa * N^gamma ────────────────────
fit <- nls(
  pangenome ~ kappa * genomes_count^gamma,
  data = rarefaction,
  start = list(kappa = max(rarefaction$pangenome[rarefaction$genomes_count == min(rarefaction$genomes_count)]), gamma = 0.5)
)

cat("\n=== Heap's Law fit ===\n")
print(summary(fit))

kappa_est <- coef(fit)["kappa"]
gamma_est <- coef(fit)["gamma"]
ci <- confint(fit)

cat("\nkappa (κ):", round(kappa_est, 2), "\n")
cat("gamma (γ, openness parameter):", round(gamma_est, 4), "\n")
cat("95% CI for gamma:", round(ci["gamma", 1], 4), "-", round(ci["gamma", 2], 4), "\n")

if (gamma_est > 0) {
  cat("\n--> gamma > 0: OPEN pangenome (still accumulating new genes with each added genome)\n")
} else {
  cat("\n--> gamma <= 0: CLOSED pangenome (approaching saturation)\n")
}

# ── Plot: data + fitted curve ──────────────────────────────────────
pred_data <- tibble(genomes_count = seq(min(rarefaction$genomes_count), max(rarefaction$genomes_count), length.out = 100)) %>%
  mutate(pangenome_pred = kappa_est * genomes_count^gamma_est)

p <- ggplot(rarefaction, aes(x = genomes_count, y = pangenome)) +
  geom_jitter(width = 0.08, alpha = 0.35, size = 1.3, color = "#B491E5") +
  geom_line(data = pred_data, aes(x = genomes_count, y = pangenome_pred), color = "grey20", linewidth = 0.8) +
  annotate("text", x = min(rarefaction$genomes_count), y = max(rarefaction$pangenome),
           label = paste0("n = ", round(kappa_est, 1), " · N^", round(gamma_est, 3)),
           hjust = 0, vjust = 1, size = 3.8, family = "Helvetica", fontface = "bold") +
  scale_x_continuous(breaks = min(rarefaction$genomes_count):max(rarefaction$genomes_count)) +
  labs(x = "Number of genomes sampled (N)", y = "Pangenome size (total gene families)",
       title = "Heap's Law fit -- pangenome openness") +
  theme_minimal(base_family = "Helvetica") +
  theme(
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    plot.title = element_text(face = "bold", size = 13)
  )

out_path <- file.path(base_dir, "heaps_law_fit.png")
ggsave(out_path, p, width = 8, height = 5.5, dpi = 300, bg = "white")
cat("\nSaved:", out_path, "\n")
