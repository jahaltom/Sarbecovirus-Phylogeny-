
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

**Script:** `ORF10_MSA.sh`  
- Extracts ORF10 from the Wuhan-Hu-1 reference genome.
- Uses `blastn` to identify homologs in all genomes.
- Extracts ORF10 hits using `bedtools` and aligns them with MAFFT.

Output: `results/ORF10_aligned.fasta`

## 3. Whole Genome Phylogeny Pipeline

**Script:** `WholeGenomeTree.sh`  
- Aligns all genomes using MAFFT.
- Filters poorly aligned regions with `trimAl`.
- Builds a bootstrapped maximum likelihood tree using `IQ-TREE` (`MFP`, SH-aLRT, and UFBoot). 1000 replicates each.
- Automatically roots the tree with `NC_019843.3` as outgroup using R/ape.

Output: `results/genomes_tree_rooted.nwk`

## 4. Tree + ORF10 Visualization

**Script:** `Plot.r`  
- Phylogenetic tree + MSA integration: Displays a rooted maximum likelihood tree alongside a multiple sequence alignment (MSA) of the ORF10 coding region.
- Clade annotation by MRCA: Automatically detects SARS-CoV-1-like and SARS-CoV-2-like clades using node-based MRCA logic and assigns branch colors.
- Support-aware and length-based node collapsing: Collapses internal nodes only if both SH-aLRT and bootstrap support values fall below user-defined thresholds (e.g., 95%), and also collapses branches with length below 0.1% of the maximum. Bootstrap values are accurately remapped to surviving nodes using bipartition (split) matching.
- Codon-level and metadata tracks: Annotates each tip with isolate metadata (host, date, location, etc.) and aligned nucleotide codons, highlighting start/stop codons and shading sequences with complete in-frame ORFs.
- Adaptive layout: Dynamically offsets the tree, metadata, and alignment panels to prevent label collision and maintain visual clarity.
- Proportional branch scale: Includes a labeled substitution scale (e.g., “0.02 substitutions/site”) with root branch adjusted to avoid artificial elongation.

Output: `ORF10_tree_MSA_group_date_codon_legend_clean_FINAL_BOOTSTRAP_FIXED.pdf`

---

## Requirements
```
mamba create -n seqAnalysis python=3.10 mafft blast bedtools trimal hyphy  iqtree entrez-direct pandas r-base r-ape -c bioconda -c conda-forge
```





