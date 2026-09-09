# Data Availability

## Newly generated data (this study)

### Philippine *L. yasudae* isolates (LyasA, LyasB)

| Isolate | Organism | BioProject | BioSample | Genome assembly (GenBank) |
|---|---|---|---|---|
| LyasA | *Leptospira yasudae* ss3a1c | PRJNA1459639 | SAMN62336287 | JCCCNU000000000 |
| LyasB | *Leptospira yasudae* ss3a1f | PRJNA1459639 | SAMN57560757 | JBYAXG000000000 |

Both isolates were deposited under the same BioProject. Raw sequencing reads are linked to their respective BioSample records; genome assembly accessions are WGS master accessions (format indicates a multi-contig draft assembly, consistent with the fragmented assembly status reported in Table 1).

## Publicly available genomes used for comparative/pangenome analysis

| Genome code | RefSeq accession | Country | Habitat |
|---|---|---|---|
| LyasC | GCF_019266585.1 | Malaysia | Soil |
| LyasD | GCF_004770815.1 | Malaysia | Water |
| LyasE | GCF_004770235.1 | Malaysia | Soil |
| LyasF | GCF_004770135.1 | Malaysia | Soil |
| LyasG | GCF_004769735.1 | Malaysia | Water |
| LyasH | GCF_004769695.1 | Malaysia | Water |
| LyasI | GCF_003545925.1 | Brazil | Soil |
| LyasJ | GCF_003545865.1 | Brazil | Soil |

All genomes retrieved from NCBI GenBank/RefSeq. Isolation country and habitat metadata compiled from associated NCBI BioSample records and original isolate-description publications (see manuscript references).

## Reference genome used for chromosome/replicon assignment

| Species | Strain | Assembly accession |
|---|---|---|
| *Leptospira interrogans* | CUDO5 | GCF_002370085.3 |

## Third-party reference databases

| Database | Version used | Source |
|---|---|---|
| Comprehensive Antibiotic Resistance Database (CARD) | latest, at time of analysis | https://card.mcmaster.ca/ |
| mobileOG-db | latest, at time of analysis | https://mobileogdb.flsi.cloud.vt.edu/ |
| *Leptospira* PubMLST scheme | latest, at time of analysis | https://pubmlst.org/organisms/leptospira-spp |
| Foldseek reference structure databases (AlphaFold/UniProt, PDB100) | latest, at time of analysis | https://search.foldseek.com/ |

> **Note:** exact database version numbers or access dates should be filled in at the time of final manuscript preparation, since these third-party resources update independently of this study's timeline.

## Source data

Numerical source data underlying the figures in the manuscript (e.g., the pairwise ANI matrix, rarefaction replicate data, per-genome/per-family mutation rates, PCA input values) are available from the corresponding author upon request, consistent with the manuscript's Data Availability statement. Where feasible, these will also be made available in this repository under `supplementary_tables/`; check there first, as this repository may be updated with additional processed data after publication.

## Code availability

All custom analysis code used to generate these outputs, including the misclustering quality-control screen, the custom dN/dS (Nei–Gojobori) implementation, and the SNP/indel/frameshift classification pipeline, is available in this repository under `scripts/`, and can be used to regenerate the same outputs from the raw/public inputs listed above. Code is released under the license specified in `LICENSE`.
