# Pruning Log — Scripts Filtered to Paper-Relevant Only

This repository was filtered from an original working set of ~102 scripts
down to 46, keeping only what supports content actually present in the
locked manuscript (4 Results subsections, Figures 1–3, Table 1, Methods,
Supplementary Methods). This log documents every exclusion decision so it
can be spot-checked — one such check already caught a real error (see
Figure 2b note below) before this was finalized.

**If anything below looks wrong, the original, unfiltered set of ~102
scripts is still available on request — nothing has been deleted, only
excluded from this filtered copy.**

## Confirmed via direct question (not inferred)

- **Figure 2b = `plot_pangenome_flower.R`.** Initially assumed this was
  superseded by later sunburst/drilldown attempts, based on the project's
  overall trajectory. This was wrong — confirmed directly that the flower
  plot is the actual final figure. `plot_pangenome_sunburst.R` and
  `plot_pangenome_drilldown.R`, along with all three Sankey attempts, are
  excluded as abandoned alternatives.

## Category-by-category reasoning

**00a/00b/00c (assembly, species/strain ID, resistance/mobile-element)** —
kept in full. Directly supports Results §1 and Table 1.

**01_coverage** — kept both mapping scripts (needed for the >1000× depth
figure in Table 1). Excluded `plot_coverage_depth.R`: no coverage-depth
figure appears in the three locked figures.

**02_ani_synteny** — kept the ANI heatmap script (Fig. 1e) and the
fastANI `--visualize` + multi-layer synteny scripts (Fig. 1f), per this
project's own established final-solution note after multiple earlier
synteny approaches were tried and abandoned. Excluded: `plot_synteny.R`,
`plot_synteny_and_ani.R`, `plot_syntenyplotter.R`,
`convert_to_syntenyplotter.py`, `build_synteny_blocks/inputs/links.py`,
`plot_whole_genome_synteny.R`, `visualize_fastani.R` (single-pair
validation only, not the final multi-genome figure), and the combined
`plot_tree_and_ani.R` / `plot_tree_and_synteny.R` (no combined panel
appears in the locked figures — Fig. 1e and 2a are separate).

**03_phylogenomics** — kept only `plot_pangenome_tree_v3.R`, inferred as
the final version given the v1→v2→v3 naming progression. **Not directly
confirmed — worth a quick check that v3 is genuinely correct.** Excluded
`plot_phylogenetic_networks.R` (NeighborNet networks never appear in the
three locked figures).

**04_pangenome_composition** — kept the flower plot (confirmed, see
above), `plot_top_gene_families_separate.R` (Fig. 2c; the non-"separate"
original version excluded as superseded), and `plot_genome_pca.R`
(Fig. 3c). Excluded all Sankey/sunburst/drilldown attempts, and
`build_cooccurrence_network.R` (gene co-occurrence network analysis never
appears in the locked figures).

**05_rarefaction_heaps_law** — kept `fit_heaps_law.R` (Fig. 2d) and
`heaps_law_convergence.R` (supports the "stabilized by N=3" claim in
Results §3/Discussion). Excluded `heaps_law_saturation.R` — this was the
extrapolation-based "~1,429 genomes" approach, which was superseded by
the convergence-based approach actually reflected in the locked text.
Excluded `plot_heaps_law.R` as redundant with `fit_heaps_law.R`'s own
built-in plotting.

**06_mutation_landscape** — kept everything supporting Fig. 3a/b/d/e/f:
SNP/indel counting, frameshift classification, family-level dN/dS, the
final length-vs-events regression script, and the validation/self-test
scripts explicitly referenced in Supplementary Methods S2. **Excluded
the per-genome dN/dS and per-genome SNP/indel scripts
(`calculate_dnds_per_genome.py`, `count_snp_indel_per_genome.py`,
`plot_dnds_per_genome.R`) — these supported an earlier three-lineage
divergence investigation (comparing genomes by sub-lineage) that does not
appear in the four locked Results subsections as currently drafted.** If
that finding should actually be part of the paper, these need to be
restored — worth double-checking this wasn't meant to be included.

**07_genome_specific_genes** — kept the identification, product-lookup,
adjacency-confirmation, and the final unified context-extraction/plotting
scripts (Fig. 2e–h). **Excluded the single-indel-event and
single-SNP-event family search scripts
(`find_single_indel_families.py`, `find_single_snp_families.py`,
`enrich_indel_families.py`, `identify_top_genes_per_metric.R`,
`check_full_contig_genes.py`) — this was a separate investigation (families
with exactly one mutational event) that does not appear in the locked
Results subsections.** Also excluded the earlier per-region extraction/plot
scripts (`extract_cluster_context.py`, `extract_lyasA_transposase_context.py`,
`extract_lyasB_remaining_context.py`, `plot_lyasA_transposase_context.R`,
`plot_lyasB_gene_cluster.R`, `plot_lyasB_remaining_genes.R`) as superseded
by the unified `extract_all_gene_contexts.py` / `plot_all_gene_contexts.R`.

**08_sequence_diagrams** — kept only the SNP-example extraction script
(Fig. 3d shows a SNP, not an indel) and the Playwright rendering
infrastructure. Excluded the earlier 2-genome pairwise illustrative
diagram scripts and the indel-alignment extraction script, as these
supported earlier exploratory figures, not the locked Fig. 3d.

**09_plasmid_investigation — excluded entirely.** Per the explicit
decision that these contigs belong to chromosome 1, not a separate
plasmid finding. No pDO5-related content appears anywhere in the locked
Results subsections.

**10_circular_genome_map** — kept in full. Figure 1a–d's description
(genome size, GC skew, GC content, contig, and feature tracks) matches
this Proksee-style circular map directly.

**11_mk_test — excluded entirely.** The McDonald-Kreitman test was never
completed and does not appear in the locked Results subsections.

**12_misc** — kept `count_gff3_features.py`; its output categories
(coding sequence, rRNA, tRNA, tmRNA, ncRNA counts) map directly onto
Table 1's annotation columns.

## Two things worth your direct confirmation

1. **`plot_pangenome_tree_v3.R` as the true final tree script** — inferred
   from version numbering, not directly confirmed the way Figure 2b was.
2. **The per-genome/three-lineage dN/dS and SNP scripts** — excluded on
   the grounds that this finding isn't in the four locked Results
   subsections. If it should be, both the scripts and the corresponding
   Results text need to be added back together.
