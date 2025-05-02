#!/bin/bash
#SBATCH --job-name=genome_phylo_bootstrap
#SBATCH --output=genome_phylo_bootstrap.out
#SBATCH --error=genome_phylo_bootstrap.err
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --mem=32G

# 1. Activate environment
echo "Setting up conda environment..."
conda activate seqAnalysis
module load R

# 2. Extract sequences from NCBI. NC_019843.3 is a Sarbecovirus outgroup. NC_045512.2 is Wuhan-Hu-1.
cat strains.tsv | awk '{print $1}' | grep -v "Taxa" | xargs -I{} efetch -db nucleotide -id {} -format fasta > all_genomes.fasta

# 3. Create output directory
mkdir -p results

# 4. Multiple sequence alignment with MAFFT
echo "Aligning full genomes with MAFFT..."
mafft --thread 8 --auto all_genomes.fasta > results/genomes_aligned.fasta

# 5. Trim poorly aligned regions using trimAl
echo "Trimming alignment with trimAl..."
trimal -in results/genomes_aligned.fasta -out results/genomes_aligned_trimmed.fasta -automated1

# 6. Build unrooted tree with IQ-TREE and bootstraps
echo "Building unrooted tree with bootstrapping using IQ-TREE..."
iqtree -s results/genomes_aligned_trimmed.fasta \
       -m MFP \
       -bb 1000 \
       -alrt 1000 \
       -nt AUTO \
       -pre results/genomes_iqtree

# 7. Root the tree using R
echo "Rooting the tree with outgroup..."
Rscript --vanilla - <<EOF
library(ape)
tree <- read.tree("results/genomes_iqtree.treefile")
outgroup_seqs <- c("NC_019843.3")
outgroup_seqs <- outgroup_seqs[outgroup_seqs %in% tree$tip.label]

if (length(outgroup_seqs) >= 1) {
  rooted_tree <- root(tree, outgroup = outgroup_seqs, resolve.root = TRUE)
  write.tree(rooted_tree, file = "results/genomes_tree_rooted.nwk")
  cat("✅ Rooted bootstrapped tree saved as results/genomes_tree_rooted.nwk\n")
} else {
  cat("⚠️ No outgroup sequences found in tree — skipping rooting.\n")
}
EOF

echo "✅ Whole-genome bootstrapped phylogenetic pipeline complete. Check the 'results' folder."












