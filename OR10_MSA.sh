#!/bin/bash
#SBATCH --job-name=genome_phylo
#SBATCH --output=genome_phylo.out
#SBATCH --error=genome_phylo.err
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --mem=32G


conda activate seqAnalysis

# 2. Extract sequences from NCBI.  NC_045512.2 is Wuhan-Hu-1.
cat strains.tsv | awk '{print $1}' | grep -v "Taxa" | grep -v NC_019843.3 |  xargs | sed 's/ /,/g' | xargs -I{} efetch -db nucleotide -id {} -format fasta > all_genomes.fasta



esearch -db nucleotide -query "NC_045512.2" | efetch -format fasta -seq_start 29558 -seq_stop 29674 > ORF10.fasta

# 3. Create output directories
mkdir -p  blastdb
mkdir -p  results

# 4. Make BLAST database
echo "Creating BLAST database..."
makeblastdb -in all_genomes.fasta -dbtype nucl -out blastdb/sarbeco_db

# 5. Run BLAST to extract ORF10 homologs
echo "Running BLAST..."
blastn -query ORF10.fasta -db blastdb/sarbeco_db -outfmt 6 -evalue 1e-5 -num_threads 8 -out results/orf10_hits.tsv

# 6. Generate BED file for matched ORF10 regions
echo "Generating ORF10 BED file..."
awk '{OFS="\t"; if ($9<$10) print $2, $9-1, $10; else print $2, $10-1, $9}' results/orf10_hits.tsv > results/orf10_coords.bed

# 7. Extract ORF10 regions from genomes
echo "Extracting ORF10 regions from genomes..."
bedtools getfasta -fi all_genomes.fasta -bed results/orf10_coords.bed -fo results/ORF10_extracted.fasta

# 8. Align with MAFFT
echo "Aligning sequences with MAFFT..."
mafft --thread 8 --auto results/ORF10_extracted.fasta > results/ORF10_aligned.fasta




echo "✅ Pipeline complete. Outputs are in the 'results' folder."
