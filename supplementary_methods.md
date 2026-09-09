# Supplementary Methods

## S1. Misclustering quality-control screen

Prior to any family-level comparative analysis, gene families were screened for evidence of erroneous clustering (e.g., merged paralogous groups) using a custom bimodal pairwise-identity test. For each gene family, all pairwise amino acid identity values among member sequences were computed. A family was flagged as a likely misclustering artifact if it satisfied both of the following conditions simultaneously: (i) minimum pairwise identity within the family < 90%, and (ii) maximum pairwise identity within the family ≥ 97%. This pattern, a family containing both near-identical and substantially divergent member pairs, is inconsistent with a single true orthologous family and is instead indicative of two or more distinct families having been erroneously merged during clustering.

Applying this screen to the persistent (core) and shell partitions identified 143 of 3,470 screened core families (4.1%) and 8 of 197 shell families (4.1%) as likely misclustering artifacts. Flagged families were excluded from all downstream family-level analyses, including SNP, indel, frameshift, and dN/dS calculations. Screen validity was confirmed using known positive and negative control cases identified during manual review of ranked outlier families.

## S2. Custom dN/dS (Nei–Gojobori) implementation

Non-synonymous (dN) and synonymous (dS) substitution counts were calculated using a self-contained implementation of the Nei–Gojobori (1986) counting method. A standard third-party codon-alignment library was evaluated initially but failed to produce valid output for the substantial majority of tested sequence pairs, and was therefore not used; the custom implementation used throughout this study does not depend on that library.

For each codon position compared between two sequences, the number of synonymous (Nd) and non-synonymous (Nd) sites and observed differences were tabulated according to the standard Nei–Gojobori path-counting procedure, averaging over all possible mutational paths where a codon differed at more than one position. Genome-wide and per-genome dN/dS values were calculated as the ratio of summed counts across all contributing codon comparisons (ΣNd / ΣN, ΣSd / ΣS), rather than as the mean of per-pair or per-codon ratios, to avoid instability arising from division by small denominators at individual sites. Prior to each analysis run, the implementation was validated against a small set of synthetic test cases with known expected outcomes (a purely synonymous codon pair, expected Nd = 0; and a purely non-synonymous codon pair, expected Nd > 0) to confirm correct behavior before application to genome-scale data.

## S3. SNP, frameshift, and in-frame indel event classification

Mutational events were classified directly from multiple sequence alignments (nucleotide, codon-partitioned) using the following column-wise procedure.

An alignment column was classified as a gap column if any aligned sequence contained a gap character at that position. Consecutive gap columns were collapsed into a single indel event, regardless of the number of columns spanned, since a contiguous run of gap columns represents one insertion or deletion event rather than one event per affected column. The length of each indel event (in nucleotides) was recorded, and events were classified as frameshift if their length was not a multiple of three, or in-frame if their length was a multiple of three.

A column was classified as a candidate single-nucleotide polymorphism (SNP) site only if it contained no gap character in any sequence (i.e., a clean, alignment-unambiguous column) and the nucleotide state was not identical across all compared sequences. Unlike indel events, each qualifying SNP column was counted as an independent event; no collapsing of adjacent variable columns was applied, since each nucleotide substitution constitutes an independent mutational event.

## S4. Cross-file gene and family identifier reconciliation

Pangenome analysis outputs referenced gene families using two distinct identifier conventions across different output files: the gene presence–absence matrix used one representative gene per family as its identifier, while per-family multiple sequence alignment files used a different representative gene, selected independently, to name the corresponding output file. Because both conventions selected a real member gene of the same underlying family, but not necessarily the same one, direct identifier matching between these files was not reliable.

To reconcile identifiers across files, a gene-to-family mapping was constructed from a third output file containing the complete membership of every individual gene copy to its family, using the identifier convention consistent with the per-family alignment files. Because the presence–absence matrix's family identifier is itself a real individual gene, it was possible to resolve the correct, alignment-file-matching family identifier by looking up that gene within the gene-to-family mapping. This reconciliation step was applied prior to all analyses requiring cross-referencing between the presence–absence matrix and per-family alignment data.

## S5. Structural homology search parameters

Candidate proteins lacking informative sequence-based functional annotation were subjected to structure-based homology search. Protein structures were predicted directly using Foldseek's integrated structure prediction and search functionality, and resulting structural matches were queried against reference structural databases. Per-residue structural prediction confidence (predicted local distance difference test, pLDDT) was used to assess the reliability of each predicted structure prior to interpreting its search results, using the following bands: pLDDT ≥ 90, very high confidence; 70–90, confident; 50–70, low confidence; < 50, very low confidence, potentially indicating a disordered region rather than a modeling failure.

Structural matches were considered high-confidence when the reported E-value was below 0.001, consistent with conventional significance thresholds for structural and sequence homology search. Matches with E-values approaching or exceeding 1 were not interpreted as evidence of true homology, as this range is consistent with matches expected by chance alone.

## S6. Heap's law model fitting and stability assessment

Pangenome openness was assessed by fitting Heap's law, expressed as pangenome size = κ·N^γ, where N is the number of genomes sampled, using nonlinear least-squares regression. Initial parameter values for optimization were set as follows: the starting value for κ was set to the observed pangenome size at the smallest sampled value of N; the starting value for γ was set to 0.5. Ninety-five percent confidence intervals for γ were obtained via profile likelihood.

To assess whether the resulting γ estimate was sensitive to the specific number of genomes included in the fitting data, the model was iteratively refit using data restricted to progressively larger subsets of the full sampled range (from the smallest subset supporting a two-parameter fit, up to the complete ten-genome dataset), and the resulting γ estimates and confidence intervals were compared across these nested subsets.

## S7. Pangenome rarefaction parameters

Pangenome and core genome rarefaction curves were generated by randomly subsampling genomes without replacement at every possible sample size from 1 to the total number of genomes analyzed, with 30 independent random permutations performed at each sample size to characterize sampling variability. All individual replicate observations, rather than per-sample-size summary statistics alone, were retained and used directly as input to the Heap's law fitting procedure described above (S6).

## S8. Detailed software commands and parameters

Commands are listed in pipeline order. Tool versions correspond to the latest stable release available at the time each analysis step was performed (see `environment.yml` in the companion code repository for the full dependency list).

**Read quality control and trimming**
```
trimmomatic PE -threads [N] -Xmx[N]g \
  input_R1.fastq.gz input_R2.fastq.gz \
  output_R1_paired.fastq.gz output_R1_unpaired.fastq.gz \
  output_R2_paired.fastq.gz output_R2_unpaired.fastq.gz \
  LEADING:30 TRAILING:30 SLIDINGWINDOW:4:30 MINLEN:50
fastqc output_R1_paired.fastq.gz output_R2_paired.fastq.gz
```
Trimming enforced a minimum Phred quality threshold of Q30 (per main Methods); explicit Java heap allocation (`-Xmx`) was required to accommodate available system memory during large paired-end trimming runs.

**Genome assembly and quality assessment**
```
spades.py --isolate -1 R1_paired.fastq.gz -2 R2_paired.fastq.gz -o assembly_out -m [N]
quast.py assembly_out/contigs.fasta -o quast_out
```
Isolate mode was used, appropriate for high-coverage, low-heterogeneity bacterial single-isolate sequencing data.

**Contamination and taxonomic screening**
```
blastn -query contigs.fasta -db nt -outfmt 6 -max_target_seqs 5 -evalue 1e-10
```
Used to screen assembled contigs for non-target (contaminant) sequence prior to downstream analysis; contigs with top hits to non-*Leptospira* taxa were flagged for removal.

**Reference-guided scaffolding**
```
ragtag.py scaffold reference.fasta contigs.fasta -o ragtag_out
```

**Genome annotation**
```
bakta --db [database_path] --output bakta_out assembly.fasta
```
Default annotation parameters were used throughout.

**Average nucleotide identity**
```
# Species/strain identity (single query vs. single reference)
fastANI -q query.fasta -r reference.fasta -o output.ani

# All-vs-all matrix
fastANI --ql genome_list.txt --rl genome_list.txt -o all_vs_all.ani

# Synteny/reciprocal mapping visualization
fastANI -q query.fasta -r reference.fasta --visualize -o output
```
Default k-mer size and fragment length parameters were used in all cases. All-vs-all ANI values were symmetrized by averaging both directions of each unordered genome pair, since fastANI's query-vs-reference calculation is not perfectly symmetric by construction.

**Whole-genome and gene-level alignment**
```
minimap2 -x asm5 reference.fasta query.fasta > output.paf     # closely related genome/chromosome assignment
minimap2 -x asm10 reference.fasta query.fasta > output.paf    # cross-species ortholog/homology search
minimap2 -x asm5 -c reference.fasta query.fasta > output.paf  # alignments requiring CIGAR strings
```
Preset stringency was selected according to expected sequence divergence: `asm5` (~0.1–5% divergence) for within-species chromosome/replicon assignment against a close reference; `asm10` (~5–10% divergence) for more divergent cross-species comparisons.

**Read mapping and coverage depth**
```
bwa index reference.fasta
bwa mem -t [N] reference.fasta R1_paired.fastq.gz R2_paired.fastq.gz | \
  samtools sort -o sorted.bam
samtools index sorted.bam
mosdepth --by 1000 sample_prefix sorted.bam
```
Coverage depth was calculated in 1 kb non-overlapping windows across the genome.

**Pangenome construction**
```
ppanggolin all --anno organisms.tsv -o pangenome_results --cpu [N]
```
Default clustering and partitioning parameters were used (MMseqs2-based clustering, statistical partitioning into persistent/shell/cloud categories).

**Pangenome rarefaction**
```
ppanggolin rarefaction -p pangenome.h5 --min 1 --max 10 --depth 30 -o rarefaction_out
```

**Core-genome phylogenomics**
```
ppanggolin msa -p pangenome.h5 --partition persistent --phylo --source dna -o msa_out
fasttree -gtr -nt persistent_genome_alignment.aln > core_genome.tree
```

**Antibiotic resistance gene screening**
```
rgi main --input_sequence contigs.fasta --output_file rgi_out --local -a DIAMOND --clean
```
The strict cut-off classification (as reported by RGI) was used to define resistance gene calls.

**Mobile genetic element screening**
```
diamond blastp --db mobileOG-db.dmnd --query proteins.faa --outfmt 6 [standard columns] \
  --evalue 1e-10 --out mobileog_output.csv
```
Screening followed the mobileOG-db recommended DIAMOND-based alignment workflow against the mobileOG-db reference database.

**Multilocus sequence typing**
```
mlst --scheme leptospira assembly.fasta
```
Allele calls were queried against the current *Leptospira* PubMLST scheme.

**Structural homology search**
Performed via Foldseek's integrated web-based structure prediction and search interface (`search.foldseek.com`), querying candidate protein sequences directly against the AlphaFold/UniProt and PDB100 reference structure databases.

**Coverage and sequence visualization**
```
# Rendering of HTML/CSS sequence diagrams
playwright screenshot --device-scale-factor=3 input.html output.png
```
A device scale factor of 3 was used throughout to ensure print-resolution (≈288 DPI) output for all rendered sequence-level illustrative figures.

**Statistical analysis and plotting**
All statistical tests, regression models, and figures were generated in R using the packages listed in `environment.yml`. Genome-level color palettes were generated deterministically using `colorspace::qualitative_hcl()` (palette = "Dark 3" or "Set 2", as specified per figure) to ensure reproducible, colorblind-conscious color assignment across all multi-genome figures.

