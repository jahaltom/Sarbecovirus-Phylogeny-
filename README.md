
# ORF10 Sarbecovirus Phylogenetic Analysis

This project reconstructs a high-resolution phylogenetic tree of Sarbecoviruses using whole-genome data and the ORF10 gene region. The pipeline incorporates automated data retrieval, alignment, filtering, tree construction, and customized visualization.

## Overview

The analysis consists of three main phases:

1. **Metadata Curation & Sequence Retrieval**  
2. **Multiple Sequence Alignment & Phylogenetic Tree Construction**  
3. **Tree Visualization with MSA and Codon Annotations**

---

## 1. Metadata Curation

**Script:** `NCBI Search.sh`  
- Queries NCBI for Sarbecovirus complete genomes (excluding SARS-CoV-2).
- Extracts isolate metadata: Accession, Isolate, Organism, Host, Date, and Location.
- Filters records with missing date/host fields and removes duplicates by `Organism`, `Host`, `Date`, and `Location`.

```bash
esearch -db nucleotide -query 'txid2509511[Organism:exp] AND "complete genome"[Title] NOT "Severe acute respiratory syndrome coronavirus 2"[Organism]' | ...
```
### Small manual edits;

- Viet Nam: Son La province, Moc Chau district, Long Sap commune, Pha nhe hamlet - cave 2 -> Viet Nam
- Add in FJ882963.1, GQ153545.1, NC_004718.3, NC_045512.2, and NC_019843.3
- Add in Group column ( "SARS-CoV", "SARS-CoV-2" "Other Sarbecovirus strains"  "Sarbecovirus Outgroup (MERS)" )

output: strains.tsv
## 2. ORF10 MSA Pipeline

**Script:** `OR10_MSA.py`  
- Extracts ORF10 from the Wuhan-Hu-1 reference genome.
- Uses `blastn` to identify homologs in all genomes.
- Extracts ORF10 hits using `bedtools` and aligns them with MAFFT.

Output: `results/ORF10_aligned.fasta`

## 3. Whole Genome Phylogeny Pipeline

**Script:** `WholeGenomeTree.py`  
- Aligns all genomes using MAFFT.
- Filters poorly aligned regions with `trimAl`.
- Builds a bootstrapped maximum likelihood tree using `IQ-TREE` (`MFP`, SH-aLRT, and UFBoot).
- Automatically roots the tree with `NC_019843.3` as outgroup using R/ape.

Output: `results/genomes_tree_rooted.nwk`

## 4. Tree + ORF10 Visualization

**Script:** `Plot.py`  
- Visualizes the rooted tree alongside the ORF10 MSA.
- Annotates clades as SARS-CoV-1-like and SARS-CoV-2-like (Colors branches) using node-based MRCA detection.
- Overlays bootstrap support.
- Displays aligned codons, isolate metadata, and highlights terminal stop codons.
- Dynamically offsets metadata and alignment positions.

Output: `ORF10_tree_MSA_group_date_codon_legend_clean_FINAL_BOOTSTRAP_FIXED.pdf`

---

## Requirements
```
mamba create -n seqAnalysis python=3.10 mafft blast bedtools trimal hyphy  iqtree entrez-direct pandas r-base r-ape -c bioconda -c conda-forge
```





