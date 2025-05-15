# Load required libraries
libs <- c("ggtree", "ggplot2", "treeio", "Biostrings", "reshape2", "ape", "ggnewscale")
for (pkg in libs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# Load tree and alignment
tree <- read.tree("genomes_tree_rooted.nwk")
aligned <- readDNAStringSet("ORF10_aligned.fasta")

# Extract SH-aLRT and Bootstrap values from node labels
sh_alrt_support <- as.numeric(sub("/.*", "", tree$node.label))
bootstrap_support <- as.numeric(sub(".*/", "", tree$node.label))

# Store node numbers before collapsing
original_node_nums <- (length(tree$tip.label) + 1):(length(tree$tip.label) + length(tree$node.label))
node_map <- data.frame(node = original_node_nums,
                       bootstrap = bootstrap_support,
                       shalrt = sh_alrt_support,
                       stringsAsFactors = FALSE)

# Collapse weak nodes if both supports are valid and below thresholds
collapse_threshold_bootstrap <- 95
collapse_threshold_shalrt <- 95
for (i in seq_len(nrow(node_map))) {
  if (!is.na(node_map$bootstrap[i]) && !is.na(node_map$shalrt[i]) &&
      node_map$bootstrap[i] < collapse_threshold_bootstrap &&
      node_map$shalrt[i] < collapse_threshold_shalrt) {
    edge_idx <- which(tree$edge[,1] == node_map$node[i])
    if (length(edge_idx) > 0) {
      tree$edge.length[edge_idx] <- 0
    }
  }
}

# Collapse into multifurcations
tree <- di2multi(tree, tol = 1e-8)

# Collapse short branches
max_len <- max(tree$edge.length)
short_branch_threshold <- max_len * 0.001
tree$edge.length[tree$edge.length < short_branch_threshold] <- 0

# Remap node labels using edge splits
require(phangorn)
original_splits <- prop.part(read.tree("genomes_tree_rooted.nwk"))
current_splits <- prop.part(tree)
matched_labels <- rep(NA, length(current_splits))
for (i in seq_along(current_splits)) {
  for (j in seq_along(original_splits)) {
    if (setequal(current_splits[[i]], original_splits[[j]])) {
      matched_labels[i] <- bootstrap_support[j]
      break
    }
  }
}
tree$node.label <- matched_labels

# Clean FASTA headers
names(aligned) <- sub(" .*", "", names(aligned))
names(aligned) <- sub(":.*", "", names(aligned))

# Add dummy outgroup
outgroup_name <- "NC_019843.3"
dummy_seq <- paste(rep("-", width(aligned)[1]), collapse = "")
aligned <- c(aligned, DNAStringSet(dummy_seq))
names(aligned)[length(aligned)] <- outgroup_name

# Match and order tree and alignment
matched_taxa <- intersect(tree$tip.label, names(aligned))
tree <- keep.tip(tree, matched_taxa)
aligned <- aligned[matched_taxa][match(tree$tip.label, names(aligned))]

# Ladderize and rescale tree
tree <- ladderize(tree, right = FALSE)
tree$edge.length <- tree$edge.length * 100
root_node <- length(tree$tip.label) + 1
tree$edge.length[which(tree$edge[,1] == root_node)] <- tree$edge.length[which(tree$edge[,1] == root_node)] / 50

# Convert to treedata and initialize plot
tree <- as.treedata(tree)
p_tree <- ggtree(tree, size = 1)

# Group clades
node_sarscov1 <- MRCA(tree, c("NC_004718.3", "OK017831.1"))
node_sarscov2 <- MRCA(tree, c("NC_045512.2", "OL674081.1"))
p_tree <- groupClade(p_tree, .node = node_sarscov2, group_name = "SARS-CoV-2-like")
p_tree <- groupClade(p_tree, .node = node_sarscov1, group_name = "SARS-CoV-1-like")

# Assign colors
p_tree$data$group <- NA
p_tree$data$group[p_tree$data$`SARS-CoV-2-like` == 1] <- "SARS-CoV-2-like"
p_tree$data$group[p_tree$data$`SARS-CoV-1-like` == 1] <- "SARS-CoV-1-like"
p_tree$data$group[is.na(p_tree$data$group)] <- "default"
p_tree$data$group_color <- ifelse(
  p_tree$data$group == "SARS-CoV-2-like", "#B8860B",
  ifelse(p_tree$data$group == "SARS-CoV-1-like", "#228B22", "black")
)
p_tree <- p_tree + aes(color = group_color) + scale_color_identity()

# Add support labels only to non-collapsed branches
p_tree$data$display_label <- NA
internal_nodes <- which(!p_tree$data$isTip)
for (i in internal_nodes) {
  label_value <- p_tree$data$label[i]
  edge_length <- p_tree$data$branch.length[i]
  if (!is.na(label_value) && !is.na(edge_length) && edge_length > 0) {
    p_tree$data$display_label[i] <- label_value
  }
}

# Load metadata
group_df <- read.table("/scr1/users/haltomj/ORF10_MSA/strains.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
group_df$Group <- trimws(group_df$Group)

# Tip metadata
tree_data <- p_tree$data
tip_data <- tree_data[tree_data$isTip, c("label", "y")]
colnames(tip_data)[1] <- "Taxa"
tip_y <- merge(tip_data, group_df, by = "Taxa", all.x = TRUE, sort = FALSE)
tip_y$label <- tip_y$Taxa

# MSA to long format
mat <- as.matrix(aligned)
df <- as.data.frame(t(mat))
df$Position <- 1:nrow(df)
long_df <- reshape2::melt(df, id.vars = "Position", variable.name = "Taxa", value.name = "Base")
long_df <- merge(long_df, tip_y, by = "Taxa")




# Codon detection
stop_codons <- c("TAA", "TAG", "TGA")
start_codons <- c("ATG")
stop_df <- data.frame()
start_df <- data.frame()
highlight_taxa <- c()

for (seqname in names(aligned)) {
  seq <- as.character(aligned[[seqname]])
  ungapped <- gsub("-", "", seq)
  if (nchar(ungapped) < 3) next
  codons <- substring(ungapped, seq(1, nchar(ungapped) - 2, 3), seq(3, nchar(ungapped), 3))
  if (sum(codons %in% stop_codons) == 1 && tail(codons, 1) %in% stop_codons) {
    highlight_taxa <- c(highlight_taxa, seqname)
  }
  pos_map <- which(strsplit(seq, "")[[1]] != "-")
  startCount <- 0
  for (i in seq_along(codons)) {
    codon <- codons[i]
    if ((i - 1) * 3 + 3 <= length(pos_map)) {
      align_pos <- pos_map[((i - 1) * 3 + 1):((i - 1) * 3 + 3)]
      if (codon %in% stop_codons)
        stop_df <- rbind(stop_df, data.frame(Taxa = seqname, Position = align_pos, CodonType = "Stop Codon"))
      if (codon %in% start_codons && startCount == 0) {
        start_df <- rbind(start_df, data.frame(Taxa = seqname, Position = align_pos, CodonType = "Start Codon"))
        startCount <- 1
      }
    }
  }
}

# Merge codons and metadata
codon_df <- rbind(stop_df, start_df)
codon_df <- merge(codon_df, tip_y[, c("Taxa", "y")], by = "Taxa")





# Set offsets
tree_max_x <- max(p_tree$data$x, na.rm = TRUE) + 8
group_bar_offset <- tree_max_x
isolate_bar_offset <- group_bar_offset + 5
host_bar_offset <- isolate_bar_offset + 8
date_bar_offset <- host_bar_offset + 7
location_bar_offset <- date_bar_offset + 7
alignment_start <- location_bar_offset + 7
long_df$Position <- long_df$Position + alignment_start
codon_df$Position <- codon_df$Position + alignment_start








# Final plot 
p_tree <- p_tree +
  geom_text2(aes(subset = (!isTip) & (!is.na(display_label)), label = display_label,
                 x = x - 1.1 , y = y + 0.4), hjust = 0, size = 2.5) +
  geom_tiplab(size = 3) +
  geom_tile(data = tip_y, aes(x = group_bar_offset, y = y, fill = Group),
            width = 1, height = 1, inherit.aes = FALSE) +
  scale_fill_manual(name = "Group", values = c(
    "SARS-CoV" = "green",
    "SARS-CoV-2" = "#f1c232",
    "Other Sarbecovirus strains" = "grey",
    "Sarbecovirus Outgroup (MERS)" = "#8c564b"
  )) +
  
  
  
  geom_rect(data = subset(tip_y, Taxa %in% highlight_taxa), aes(ymin = y - 0.5, ymax = y + 0.5, xmin = 0, xmax = max(long_df$Position) + 5), fill = "yellow", alpha = 0.3, inherit.aes = FALSE) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = codon_df, aes(x = Position, y = y, fill = CodonType), width = 1, height = 1, inherit.aes = FALSE) +
  scale_fill_manual(name = "Codon Type", values = c("Stop Codon" = "#fcaeae", "Start Codon" = "#4169e1")) +
  geom_text(data = tip_y, aes(x = isolate_bar_offset, y = y, label = Isolate), size = 2.5, inherit.aes = FALSE) +
  geom_text(data = tip_y, aes(x = host_bar_offset, y = y, label = Host), size = 2.5, inherit.aes = FALSE) +
  geom_text(data = tip_y, aes(x = date_bar_offset, y = y, label = Date), size = 2.5, inherit.aes = FALSE) +
  geom_text(data = tip_y, aes(x = location_bar_offset, y = y, label = Location), size = 2.5, inherit.aes = FALSE) +
  ggnewscale::new_scale_color() +
  geom_text(data = long_df, aes(x = Position, y = y, label = Base, color = "black"), size = 3, family = "mono", inherit.aes = FALSE) +
  xlim(0, max(long_df$Position) + 5) +
  
  geom_treescale(x = 5, y = 100, width = 2.0,
                 label = "0.02 substitutions/site", linesize = 1) +
  scale_color_identity() +
  
  
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.position = "right"
  )

# Save output
ggsave("ORF10_tree_MSA.png",
       plot = p_tree, width = 49, height = 32)
cat("✅ Tree figure saved\n")

