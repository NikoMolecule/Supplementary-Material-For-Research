setwd("D:/AlphaAmylase/Data")
FAM163A <- readLines("FAM163A_mammals.fasta")

library(msa)
library(ape)
library(ggtree)
library(ggplot2)
library(phangorn)
library(gridExtra)
library(bio3d)

plot_tree <- function(tree_plot, title_plot, max_x) {
  g <- ggtree(tree_plot, colour = "#444444", size = 0.6, ladderize = TRUE) +
    
    geom_tiplab(size = 3, colour = "black", align = TRUE, linesize = 0.3) +
    
    geom_nodepoint(size = 2, colour = "#c7254e", alpha = 0.8) +
    
    theme_tree2() + 
    
    xlim(0, max_x) +
    
    labs(title = title_plot) +
    
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 10, unit = "pt"),
      legend.position = "none"
    )
  
  return(g)
}

# ------------------ GETTING DATA FROM FASTA FILE ---------------------------
namesList = c()
sequences = c()

sequence = ""
for (line in FAM163A) {
  if(substring(line, 1, 1) == ">") {
    if(sequence != ""){
      sequences <- append(sequences, sequence)
      sequence = ""
    }
    name = ""
    t = FALSE
    for (i in 2:nchar(line)) {
      if(substring(line, i - 1, i - 1) == "=") {
        t = TRUE
      }
      if(substring(line, i, i) == "O" && substring(line, i + 1, i + 1) == "X" && substring(line, i + 2, i + 2) == "="){
        t = FALSE
        break
      }
      if(t == TRUE) {
        name = paste0(name, substring(line, i, i))
      }
    }
    
    namesList <- append(namesList, name)
  } else {
    sequence = paste0(sequence, line)
  }
}

if(sequence != ""){
  sequences <- append(sequences, sequence)
  sequence = ""
}

names(sequences) = namesList


# ------------ MSA USING MUSCLE ---------------------

myAAseq <- AAStringSet(sequences)

myAlignment <- msa(myAAseq, "Muscle", type = "protein")

print(myAlignment)

msaPrettyPrint(myAlignment, output="tex", showLogo="top", consensusColors = "Gray",
               logoColors="rasmol", shadingMode="similar", showConsensus="bottom",
               showLegend=TRUE, askForOverwrite=FALSE)


aln_matrix <- as.matrix(myAlignment)

# ------------- GENERATING CONSENSUS PLOT -----------------------------
raw_matrix <- aln_matrix

gap_fraction <- colMeans(raw_matrix == "-")

keep_cols <- gap_fraction < 0.5

clean_matrix <- raw_matrix[, keep_cols]

aln_matrix <- clean_matrix

last_row_index <- nrow(aln_matrix)
consensus_seq_vector <- aln_matrix[last_row_index, ]

consensus_string <- paste(consensus_seq_vector, collapse = "")

write(consensus_string, file = "consensus_sequence.txt")
print("Success! Consensus sequence saved to 'consensus_sequence.txt'")

real_species <- aln_matrix[-last_row_index, ] 
scores <- conserv(real_species, method="identity")
print(scores)

plot_data <- data.frame(
  Position = 1:length(scores),
  Score = scores
)

ggplot(plot_data, aes(x = Position, y = Score)) +
  geom_area(fill = "dodgerblue", alpha = 0.6) +
  geom_line(color = "blue", size = 0.5) +
  theme_minimal() +
  labs(title = "Conservation Landscape of FAM163A",
       y = "Conservation Score (Identity)",
       x = "Amino Acid Position") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 12))

write.csv(plot_data, "FAM163A_conservation_scores.csv", row.names = FALSE)


# ------------ DISTANCE MATRIX: FITCH SIMILLARITY MATRIX -----------------
phyDat_alignment <- as.phyDat(aln_matrix, type = "AA")

trimmed_matrix_check <- as.character(phyDat_alignment)
trimmed_identity <- seqidentity(as.matrix(trimmed_matrix_check))

pdistance <- dist.hamming(phyDat_alignment, ratio = TRUE)
mldistance <- dist.ml(phyDat_alignment, model = "LG")
# ------------ NEIGHBOUR JOINING PHYLOGENY -------------------------------
bioNJ_tree <- bionj(mldistance)
bioNJ_tree$edge.length[which(bioNJ_tree$edge.length<  0)] = 0
bioNJ_tree = midpoint(multi2di(bioNJ_tree))


plot_tree(bioNJ_tree, "FAM163A Phylogeny in Mammals", 1)
# ------------------- BOOTSTRAPPING -------------------------

bs_nj <- bootstrap.phyDat(phyDat_alignment, FUN = function(x) bionj(dist.ml(x, model = "LG")), bs = 1000)

bs_confidence_nj <- addConfidences(bioNJ_tree, bs_nj)
# ----------------- PLOTTING BOOTSTRAP VALUED TREE -----------------
na_index <- which(is.na(bs_confidence_nj$node.label))
bs_confidence_nj$node.label[na_index] <- 1


color_palette <- colorRampPalette(c("red3", "mediumseagreen"))(100)


color_indices <- round(bs_confidence_nj$node.label * 100)
node_colors <- color_palette[color_indices]


tip_colors <- rep("black", length(bs_confidence_nj$tip.label))

# Plot the tree with bootstrap node support
g <- ggtree(bs_confidence_nj, layout = "dendrogram", ladderize = TRUE)
g <- g + geom_tiplab(size = 4, color = tip_colors) 
g <- g + theme(plot.margin = margin(t = 0, r = 0, b = 150, l = 0, unit = "pt"), legend.position = "top")
g <- g + geom_nodepoint(size = 4, color = node_colors)
g <- g + labs(title = "Bootstrapped FAM163A Phylogenetic Tree")
g <- g + theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))


color_data = data.frame(x = 1:100, y = 1)
color_data$colors = color_palette

# Create a color bar plot using ggplot2
g1 = ggplot(color_data, aes(x = x, y = y, fill = colors)) +
  geom_bar(stat = "identity", width = 1, position = "identity") + scale_fill_manual(values = rev(color_data$colors)) +
  labs(title = NULL,
       x = "% of node support",
       y = NULL) +
  theme_minimal() + theme(axis.text.y = element_blank(),
                          axis.title.y = element_blank(),
                          panel.grid = element_blank()) 
g1 = g1 + guides(fill = "none")



grid.arrange(g, g1, nrow = 2, ncol = 1, heights = c(0.8, 0.2))
