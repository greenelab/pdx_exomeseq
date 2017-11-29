# Gregory Way - 2017
# pdx_exomeseq
# scripts/viz/oncoprint_similarity.R
#
# Output:
# Oncoprint diagrams and similarity heatmaps for all replicates and consensus
#
# Data processing and specific plotting functions modified from vignette:
# https://bioconductor.org/packages/3.7/bioc/vignettes/ComplexHeatmap/inst/doc/s8.oncoprint.html

library(gplots)
library(RColorBrewer)
library(ComplexHeatmap)

# Input Files
replicate_oncoprint_file <- file.path("results", "oncoprint_replicates.tsv")
consensus_oncoprint_file <- file.path("results", "oncoprint_consensus.tsv")

replicate_sim_file <- file.path("results", "cosmic_similarity_replicates.tsv")
consensus_sim_file <- file.path("results", "cosmic_similarity_consensus.tsv")

# Output Files
replicate_oncoprint_out <- file.path("figures", "oncoprint_replicates.pdf")
consensus_oncoprint_out <- file.path("figures", "oncoprint_consensus.pdf")

replicate_sim_out <- file.path("figures", "cosmic_similarity_replicates.pdf")
consensus_sim_out <- file.path("figures", "cosmic_similarity_consensus.pdf")

# Define Constants and Functions
col = c("MUT" = "#008000", "AMP" = "red", "HOMDEL" = "blue")

process_oncoprint <- function(onco_file) {
  mat <- read.table(onco_file, header = TRUE, stringsAsFactors = FALSE,
                    sep = "\t")
  mat[is.na(mat)] = ""
  rownames(mat) = mat[, 1]
  mat = mat[, -1]
  mat = mat[, -ncol(mat)]
  mat = t(as.matrix(mat))
  return(mat)
}

process_similarity <- function(sim_file) {
  cosmic_sim_df <- readr::read_tsv(sim_file)
  cosmic_dist_df <- data.matrix(dist(cosmic_sim_df, diag = TRUE, upper = TRUE))
  cosmic_dist_df <- data.matrix(cor(t(cosmic_sim_df[, 2:ncol(cosmic_sim_df)])))
  rownames(cosmic_dist_df) <- colnames(cosmic_dist_df) <- cosmic_sim_df$Case.ID
  return(cosmic_dist_df)
}

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  HOMDEL = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
              gp = gpar(fill = "blue", col = NA))
  },
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h - unit(0.5, "mm"),
              gp = gpar(fill = "red", col = NA))
  },
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h*0.33,
              gp = gpar(fill = "#008000", col = NA))
  }
)

# Process input oncoprint matrices
mat_rep <- process_oncoprint(replicate_oncoprint_file)
mat_con <- process_oncoprint(consensus_oncoprint_file)
 
# Output OncoPrint Visualizations
pdf(replicate_oncoprint_out, height = 10, width = 19)
oncoPrint(mat_rep,
          get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun,
          col = col,
          show_column_names = TRUE,
          column_title = "PDX WES Mutations (Replicates)",
          axis_gp = gpar(fontsize = 4),
          heatmap_legend_param = list(title = "Alternations",
                                      at = c("AMP", "HOMDEL", "MUT"), 
                                      labels = c("Amplification",
                                                 "Deep deletion", "Mutation")))
dev.off()

pdf(consensus_oncoprint_out, height = 10, width = 14)
oncoPrint(mat_con,
          get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun,
          col = col,
          show_column_names = TRUE,
          column_title = "PDX WES Mutations (Consensus)",
          axis_gp = gpar(fontsize = 4),
          heatmap_legend_param = list(title = "Alternations",
                                      at = c("AMP", "HOMDEL", "MUT"), 
                                      labels = c("Amplification",
                                                 "Deep deletion", "Mutation")))
dev.off()

# Process similarity matrix data
sim_rep <- process_similarity(replicate_sim_file)
sim_con <- process_similarity(consensus_sim_file)

# Output Similarity Matrix Visualizations
pdf(replicate_sim_out, height = 8, width = 8)
heatmap.2(sim_rep,
          labCol = FALSE,
          trace = "none",
          revC = TRUE,
          col = colorRampPalette(c("red", "black", "green"))(n = 200),
          main = "COSMIC Profile Similarity",
          cexRow = 0.4)
dev.off()

pdf(consensus_sim_out, height = 8, width = 8)
heatmap.2(sim_con,
          labCol = FALSE,
          trace = "none",
          revC = TRUE,
          col = colorRampPalette(c("red", "black", "green"))(n = 200),
          main = "COSMIC Profile Similarity",
          cexRow = 0.8)
dev.off()
