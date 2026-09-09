# ============================================================
# At what sample size does the gamma ESTIMATE itself stabilize below
# 0.3? Refits Heap's Law using only the rarefaction data through each
# possible max N (2 through 10), showing how the estimate (and its
# uncertainty) evolves as more genomes are added -- directly answering
# "would N genomes have been enough to conclude gamma < 0.3" using
# data we already have, no extrapolation required.
# ============================================================

library(tidyverse)

base_dir <- "~/Downloads/macrogen_1c/pangenome/pangenome_results/rarefaction"
rarefaction <- read_csv(file.path(base_dir, "rarefaction.csv"), show_col_types = FALSE)

results <- tibble()

for (max_n in 2:max(rarefaction$genomes_count)) {
  subset_data <- rarefaction %>% filter(genomes_count <= max_n)

  fit <- tryCatch(
    nls(pangenome ~ kappa * genomes_count^gamma, data = subset_data,
        start = list(kappa = max(subset_data$pangenome[subset_data$genomes_count == min(subset_data$genomes_count)]), gamma = 0.5)),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    results <- bind_rows(results, tibble(max_n = max_n, gamma = NA, ci_low = NA, ci_high = NA, below_0.3 = NA))
    next
  }

  gamma_est <- coef(fit)["gamma"]
  ci <- tryCatch(confint(fit), error = function(e) matrix(NA, 2, 2, dimnames = list(c("kappa","gamma"), c("2.5%","97.5%"))))

  results <- bind_rows(results, tibble(
    max_n = max_n,
    n_rows = nrow(subset_data),
    gamma = gamma_est,
    ci_low = ci["gamma", 1],
    ci_high = ci["gamma", 2],
    below_0.3 = gamma_est < 0.3,
    ci_entirely_below_0.3 = !is.na(ci["gamma", 2]) && ci["gamma", 2] < 0.3
  ))
}

cat("=== Gamma estimate as a function of how many genomes' worth of data went into the fit ===\n\n")
print(results, n = Inf)

cat("\n=== Interpretation ===\n")
first_stable <- results %>% filter(ci_entirely_below_0.3) %>% slice(1)
if (nrow(first_stable) > 0) {
  cat(sprintf("Smallest max_n where the estimate AND its full 95%% CI sit below 0.3: max_n = %d\n", first_stable$max_n))
} else {
  cat("No max_n in this range gives a CI entirely below 0.3 -- the estimate may fluctuate\n")
  cat("or remain uncertain until more genomes are included.\n")
}

cat("\nSpecifically for max_n = 6:\n")
row6 <- results %>% filter(max_n == 6)
print(row6)
