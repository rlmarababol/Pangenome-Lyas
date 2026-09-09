# Pangenomic analysis of *Leptospira yasudae* reveals habitat-driven core and accessory genome divergence

This repository contains the analysis code, scripts, and processed data supporting the manuscript:

> Marababol, R.J.L. & Rivera, W.L. [Year]. Pangenomic analysis of *Leptospira yasudae* reveals habitat-driven core and accessory genome divergence. *[Journal]*. [DOI once available]

## Overview

This repository documents the full computational pipeline used to:
- Confirm species and strain identity of two Philippine *L. yasudae* isolates via genome assembly, ANI, and MLST
- Construct a ten-genome pangenome (2 Philippine isolates + 8 publicly available genomes) and assess pangenome openness
- Identify genome-specific genes and characterize candidate horizontally acquired defense systems
- Quantify and compare mutational signatures (SNP, frameshift, in-frame indel, dN/dS) between core and shell gene partitions

## Repository structure

```
.
├── README.md                          # this file
├── DATA_AVAILABILITY.md               # accession numbers and data sources
├── CITATION.cff                       # machine-readable citation metadata
├── environment.yml                    # conda environment specification
├── LICENSE                            # code license
├── supplementary_methods.md           # full technical/algorithmic detail
├── scripts/
│   ├── 00a_assembly_qc/                # trimming, assembly, QC, scaffolding, annotation
│   ├── 00b_species_strain_id/          # ANI species/strain confirmation, MLST
│   ├── 00c_resistance_mobile/          # CARD resistance genes, mobileOG-db screening
│   ├── 01_coverage/                   # read mapping, mosdepth coverage analysis
│   ├── 02_ani_synteny/                # FastANI, minimap2, synteny visualization
│   ├── 03_phylogenomics/              # core-genome tree construction, phylogenetic networks
│   ├── 04_pangenome_composition/      # pangenome partitioning, composition figures, PCA
│   ├── 05_rarefaction_heaps_law/      # pangenome rarefaction and Heap's law fitting
│   ├── 06_mutation_landscape/         # SNP/frameshift/indel classification, dN/dS
│   ├── 07_genome_specific_genes/      # unique gene identification, structural homology search
│   ├── 08_sequence_diagrams/          # illustrative sequence-level alignment figures
│   ├── 09_plasmid_investigation/      # pDO5 homology and replicon assignment analysis
│   ├── 10_circular_genome_map/        # Proksee-style circular genome visualization
│   ├── 11_mk_test/                    # McDonald-Kreitman test (outgroup ortholog search)
│   └── 12_misc/                       # GFF3 feature census and other supporting scripts
└── supplementary_tables/              # supplementary data tables referenced in the manuscript
```

Scripts are organized thematically rather than in strict pipeline order, since many analyses (e.g., genome-specific gene characterization) were iterative and built on outputs from multiple earlier stages. Each script's header comment describes its specific purpose, inputs, and outputs.

## Reproducing this analysis

1. Clone this repository and create the conda environment:
   ```
   git clone [repository URL]
   cd [repository name]
   conda env create -f environment.yml
   conda activate [environment name]
   ```
2. Raw sequencing data and reference genomes must be downloaded separately; see `DATA_AVAILABILITY.md` for accession numbers and sources.
3. Scripts are numbered in the order they were used in the analysis pipeline; each subdirectory contains its own short README describing script-specific inputs and outputs.

## Data and code availability

See `DATA_AVAILABILITY.md` for full accession numbers, database versions, and third-party data sources used in this study.

## Citation

If you use this code or data, please cite the manuscript above. See `CITATION.cff` for machine-readable citation metadata.

## Contact

For questions regarding this repository, please contact [corresponding author email] or open an issue on this repository.
