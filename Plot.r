# Load required libraries
libs <- c("ggtree", "ggplot2", "treeio", "Biostrings", "reshape2", "ape", "ggnewscale")
for (pkg in libs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# Load tree and alignment
tree <- read.tree("genomes_tree_rooted.nwk")
aligned <- readDNAStringSet("ORF10_aligned.fasta")

# Extract SH-aLRT and Bootstrap values
if (!is.null(tree$node.label)) {
  sh_alrt_support <- as.numeric(sub("/.*", "", tree$node.label))
  bootstrap_support <- as.numeric(sub(".*/", "", tree$node.label))
} else {
  stop("Tree does not contain bootstrap values!")
}

# Collapse weak nodes
collapse_threshold_bootstrap <- 95
collapse_threshold_shalrt <- 95
weak_nodes <- which(bootstrap_support < collapse_threshold_bootstrap | sh_alrt_support < collapse_threshold_shalrt)

if (length(weak_nodes) > 0) {
  for (node in weak_nodes) {
    edge_idx <- which(tree$edge[,1] == node + length(tree$tip.label))
    if (length(edge_idx) > 0) {
      tree$edge.length[edge_idx] <- 0
    }
  }
}

# Collapse into multifurcations
tree <- di2multi(tree, tol = 1e-8)

# Update node labels cleanly
# Keep bootstrap values for internal tree storage
# Display separately for plot
real_node_label <- ifelse(
  (bootstrap_support >= collapse_threshold_bootstrap | sh_alrt_support >= collapse_threshold_shalrt),
  bootstrap_support,
  NA
)
tree$node.label <- real_node_label

# Create display labels for plotting
pretty_display_labels <- ifelse(
  (bootstrap_support >= collapse_threshold_bootstrap | sh_alrt_support >= collapse_threshold_shalrt),
  paste0(round(bootstrap_support)),    #round(sh_alrt_support),"%"
  NA
)

# Clean FASTA headers
names(aligned) <- sub(" .*", "", names(aligned))
names(aligned) <- sub(":.*", "", names(aligned))



# 1. Set outgroup name
outgroup_name <- "NC_019843.3"
# 2. Create a dummy sequence
dummy_width <- width(aligned)[1]    # Get the width (should be 119)
dummy_seq <- paste(rep("-", dummy_width), collapse = "")
# 3. Create a DNAStringSet object
dummy_dna <- DNAStringSet(dummy_seq)
names(dummy_dna) <- outgroup_name
# 4. Add dummy sequence to aligned
aligned <- c(aligned, dummy_dna)

# Match tree and alignment
matched_taxa <- intersect(tree$tip.label, names(aligned))
tree <- keep.tip(tree, matched_taxa)
aligned <- aligned[matched_taxa][match(tree$tip.label, names(aligned))]

# Ladderize
# Ladderize
tree <- ladderize(tree, right = FALSE)

# Scale branch lengths
tree$edge.length <- tree$edge.length * 13 + 2

# 🔧 Shorten the root branch (edge from node 0 to root)
# Find the edge leading to the root node
root_node <- length(tree$tip.label) + 1
root_edge_index <- which(tree$edge[,1] == root_node)

# Shrink its length (e.g., divide by 5)
tree$edge.length[root_edge_index] <- tree$edge.length[root_edge_index] / 5




# Convert tree to treedata
tree <- as.treedata(tree)
# === Plotting ===
p_tree <- ggtree(tree, size=1)
#Find Sars-Cov1/2 nodes
node_sarscov1 =MRCA(tree, c("NC_004718.3","OK017831.1"))
node_sarscov2 = MRCA(tree, c("NC_045512.2","OK017805.1"))
# Group clades
p_tree <- groupClade(p_tree, .node = node_sarscov2, group_name = "SARS-CoV-2-like")
p_tree <- groupClade(p_tree, .node = node_sarscov1, group_name = "SARS-CoV-1-like")
# Correct group column creation
p_tree$data$group <- NA
p_tree$data$group[p_tree$data$`SARS-CoV-2-like` == 1] <- "SARS-CoV-2-like"
p_tree$data$group[p_tree$data$`SARS-CoV-1-like` == 1] <- "SARS-CoV-1-like"
p_tree$data$group[is.na(p_tree$data$group)] <- "default"
# NEW! Real color names matching the group
p_tree$data$group_color <- ifelse(
  p_tree$data$group == "SARS-CoV-2-like", "#B8860B",
  ifelse(p_tree$data$group == "SARS-CoV-1-like", "#228B22", "black")
)
# Apply color mapping
p_tree <- p_tree + aes(color = group_color) +scale_color_identity()




# Load metadata
group_df <- read.table("/scr1/users/haltomj/ORF10_MSA/strains.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
group_df$Group <- trimws(group_df$Group)
group_df$Group[group_df$Group == "Other SARSr-CoV strains"] <- "Other SARSr-CoV strains"

# Plot tree
#p_tree <- ggtree(tree, size = 1)

# Map pretty display labels
p_tree$data$display_label <- NA
internal_nodes <- which(!p_tree$data$isTip)
p_tree$data$display_label[internal_nodes] <- pretty_display_labels

# Map bootstrap (internal)
p_tree$data$bootstrap <- NA
p_tree$data$bootstrap[internal_nodes] <- real_node_label

# Prepare tip labels
tree_data <- p_tree$data
tip_data <- tree_data[tree_data$isTip, c("label", "y")]
colnames(tip_data)[1] <- "Taxa"
tip_y <- merge(tip_data, group_df, by = "Taxa", all.x = TRUE, sort = FALSE)
tip_y$label <- tip_y$Taxa

# Prepare MSA
mat <- as.matrix(aligned)
df <- as.data.frame(t(mat))
df$Position <- 1:nrow(df)
long_df <- reshape2::melt(df, id.vars = "Position", variable.name = "Taxa", value.name = "Base")
long_df <- merge(long_df, tip_y, by = "Taxa")





# Set x-offsets
tree_max_x <- max(p_tree$data$x, na.rm = TRUE) +8
group_bar_offset <- tree_max_x 
isolate_bar_offset <- tree_max_x + 5
host_bar_offset <- tree_max_x + 13
date_bar_offset <- tree_max_x + 20
location_bar_offset <- tree_max_x + 27 
alignment_start <- tree_max_x + 30
long_df$Position <- long_df$Position + alignment_start

# Detect codons
stop_codons <- c("TAA", "TAG", "TGA")
start_codons <- c("ATG")
stop_df <- data.frame()
start_df <- data.frame()
highlight_taxa <- c()

for (seqname in names(aligned)) {
  seq <- as.character(aligned[[seqname]])
  ungapped <- gsub("-", "", seq)

  # SKIP if no ungapped sequence (e.g., dummy outgroup)
  if (nchar(ungapped) < 3) {
    next
  }

  codons <- substring(ungapped, seq(1, nchar(ungapped) - 2, 3), seq(3, nchar(ungapped), 3))
  stop_count <- sum(codons %in% stop_codons)
  last_codon <- tail(codons, 1)
  if (stop_count == 1 && last_codon %in% stop_codons) {
    highlight_taxa <- c(highlight_taxa, seqname)
  }
  pos_map <- which(strsplit(seq, "")[[1]] != "-")
  startCount <- 0
  for (i in seq_along(codons)) {
    codon <- codons[i]
    if ((i - 1) * 3 + 3 <= length(pos_map)) {
      align_pos <- pos_map[((i - 1) * 3 + 1):((i - 1) * 3 + 3)]
      if (codon %in% stop_codons) {
        stop_df <- rbind(stop_df, data.frame(Taxa = seqname, Position = align_pos, CodonType = "Stop Codon"))
      }
      if (codon %in% start_codons & startCount == 0) {
        start_df <- rbind(start_df, data.frame(Taxa = seqname, Position = align_pos, CodonType = "Start Codon"))
        startCount <- 1
      }
    }
  }
}


codon_df <- rbind(stop_df, start_df)
codon_df <- merge(codon_df, tip_y[, c("Taxa", "y")], by = "Taxa")
codon_df$Position <- codon_df$Position + alignment_start

long_df$Highlight <- ifelse(long_df$Taxa %in% highlight_taxa, TRUE, FALSE)






# Plot
p_tree <- p_tree +


  geom_text2(
    aes(subset = (!isTip) & (!is.na(display_label)),
        label = display_label,
        x = x - 1.5, y = y + 1),
    hjust = 0,
    size = 2.5
  ) +
  geom_rect(
    data = subset(tip_y, Taxa %in% highlight_taxa),
    aes(ymin = y - 0.5, ymax = y + 0.5, xmin = 0, xmax = max(long_df$Position) + 5),
    fill = "yellow",
    alpha = 0.3,
    inherit.aes = FALSE
  ) +
  geom_tiplab(size = 3) +
  geom_tile(data = tip_y, aes(x = group_bar_offset, y = y, fill = Group), width = 1, height = 1,inherit.aes = FALSE) +
  scale_fill_manual(name = "Group", values = c(
    "SARS-CoV" = "green",
    "SARS-CoV-2" = "#f1c232",
    "Other Sarbecovirus strains" = "grey",
    "Sarbecovirus Outgroup (MERS)" = "#8c564b"

    
  )) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = codon_df, aes(x = Position, y = y, fill = CodonType), width = 1, height = 1,inherit.aes = FALSE) +
  scale_fill_manual(name = "Codon Type", values = c("Stop Codon" = "#fcaeae", "Start Codon" = "#4169e1")) +
  geom_text(data = tip_y, aes(x = isolate_bar_offset, y = y, label = Isolate), size = 2.5,inherit.aes = FALSE) +
  geom_text(data = tip_y, aes(x = host_bar_offset, y = y, label = Host), size = 2.5,inherit.aes = FALSE) +
  geom_text(data = tip_y, aes(x = date_bar_offset, y = y, label = Date), size = 2.5,inherit.aes = FALSE) +
  geom_text(data = tip_y, aes(x = location_bar_offset, y = y, label = Location), size = 2.5,inherit.aes = FALSE) +
  ggnewscale::new_scale_color() +
  geom_text(
    data = long_df,
    aes(x = Position, y = y, label = Base, fontface = "bold", color = "black"),
    size = 3,
    family = "mono",inherit.aes = FALSE
  ) +
 
  scale_color_identity() +
  xlim(0, max(long_df$Position) + 5) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(20, 20, 20, 20),
    legend.position = "right",
    legend.box = "vertical"
  )

# Save output
ggsave("ORF10_tree_MSA_group_date_codon_legend_clean_FINAL_BOOTSTRAP_FIXED.pdf", plot = p_tree, width = 49, height = 32)
cat("\U00002705 PDF saved: ORF10_tree_MSA_group_date_codon_legend_clean_FINAL_BOOTSTRAP_FIXED.pdf\n")
