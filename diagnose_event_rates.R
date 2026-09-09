library(tidyverse)

screen <- read_csv("~/Downloads/macrogen_1c/pangenome/pangenome_results/family_identity_screen.csv",
                    show_col_types = FALSE)
clean_ids <- screen %>% filter(!flagged_misclustered) %>% select(category, family_id)

df <- read_csv("~/Downloads/macrogen_1c/pangenome/pangenome_results/snp_indel_events.csv",
               show_col_types = FALSE) %>%
  inner_join(clean_ids, by = c("category", "family_id"))

cat("Excluded", nrow(screen %>% filter(flagged_misclustered)), "likely-misclustered families\n\n")

cat("=== Alignment length distribution (all families) ===\n")
print(summary(df$aln_len))

cat("\n=== events_per_kb: mean vs median, by category ===\n")
df %>%
  group_by(category) %>%
  summarise(
    mean_epk   = mean(events_per_kb),
    median_epk = median(events_per_kb),
    q25        = quantile(events_per_kb, 0.25),
    q75        = quantile(events_per_kb, 0.75),
    min_epk    = min(events_per_kb),
    max_epk    = max(events_per_kb),
    .groups = "drop"
  ) %>%
  print()

cat("\n=== 10 shortest alignments and their rates (checking for tiny-denominator inflation) ===\n")
df %>%
  arrange(aln_len) %>%
  select(category, family_id, aln_len, snp, indel_events, events_per_kb) %>%
  head(10) %>%
  print(n = 10)

cat("\n=== Does short alignment length predict an inflated rate? (Spearman correlation) ===\n")
print(cor.test(df$aln_len, df$events_per_kb, method = "spearman"))

cat("\n=== What fraction of families have aln_len < 150bp? ===\n")
df %>%
  group_by(category) %>%
  summarise(pct_short = mean(aln_len < 150) * 100, .groups = "drop") %>%
  print()

cat("\n=== n_seqs (genomes present) distribution by category ===\n")
df %>%
  group_by(category) %>%
  summarise(min_n = min(n_seqs), median_n = median(n_seqs),
            mean_n = mean(n_seqs), max_n = max(n_seqs), .groups = "drop") %>%
  print()

cat("\n=== Fraction of families present in only 2 genomes (auto-zeroed by MIN_MINOR_COUNT=2) ===\n")
df %>%
  group_by(category) %>%
  summarise(pct_n2 = mean(n_seqs == 2) * 100, .groups = "drop") %>%
  print()

cat("\n=== Top 10 highest-rate Core families remaining after exclusion (sanity check) ===\n")
df %>%
  filter(category == "Core") %>%
  arrange(desc(events_per_kb)) %>%
  select(family_id, n_seqs, aln_len, snp, indel_events, events_per_kb) %>%
  head(10) %>%
  print(n = 10)
