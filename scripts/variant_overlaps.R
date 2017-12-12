# Gregory Way - 2017
# pdx_exomeseq
# scripts/viz/variant_overlaps.R
#
# Visualize overlapping called variants
#
# Output:
# Three Venn Diagrams displaying variants called across replicates (all
# variants and COSMIC variants), and COSMIC variants across sample passages.
# Also output is a summary dataframe of all common COSMIC variants called in
# each replicate.

library(dplyr)
library(gridExtra)
library(VennDiagram)

variant_dir <- file.path("results", "processed_vcfs")
variant_files <- list.files(variant_dir)

base_samples <- sapply(strsplit(variant_files, "_"), `[`, 1)

# Obtain Venn diagrams across replicates for all variants and COSMIC variants
all_variant_venns <- list()
cosmic_variant_venns <- list()
all_var_df <- list()
replicate_colors <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072")
for (samp in unique(base_samples)) {
  select_variant_files <- variant_files[grepl(samp, variant_files)]
  sample_variant_df <- c()
  for (select_file in select_variant_files) {
    variant_df <- suppressMessages(readr::read_tsv(file.path(variant_dir,
                                                             select_file)))
    variant_df$replicate <- strsplit(select_file, "_")[[1]][3]
    sample_variant_df <- rbind(sample_variant_df, variant_df)
  }
  sample_variant_df$variant_id <-
    paste(sample_variant_df$Chr, sample_variant_df$Start,
          sample_variant_df$End, sample_variant_df$Ref,
          sample_variant_df$Alt, sample_variant_df$Gene.refGene, sep = "_")
  sample_variant_df$sample_id <- samp
  sample_variant_df$base_id <- strsplit(samp, "-")[[1]][1]
  
  # Process file and output two venn diagrams
  variant_acast_full <- reshape2::acast(sample_variant_df,
                                        variant_id ~ replicate,
                                        value.var = "variant_id")
  variant_list_full <- lapply(as.data.frame(variant_acast_full), 
                              function(x) x[!is.na(x)])
  venn_plot_full_obj <- venn.diagram(x = variant_list_full,
                                     main = samp,
                                     main.pos = c(0.5, 1),
                                     main.cex = 2,
                                     fill = replicate_colors,
                                     filename = NULL)
  all_variant_venns[[samp]] <- venn_plot_full_obj
  
  # COSMIC variants
  cosmic_variant_df <- sample_variant_df %>% dplyr::filter(cosmic70 != ".")
  variant_acast_cosmic <- reshape2::acast(cosmic_variant_df,
                                          variant_id ~ replicate,
                                          value.var = "variant_id")
  variant_list_cosmic <- lapply(as.data.frame(variant_acast_cosmic),
                                function(x) x[!is.na(x)])
  venn_plot_cosmic_obj <- venn.diagram(x = variant_list_cosmic,
                                       main = samp,
                                       main.pos = c(0.5, 1),
                                       main.cex = 2,
                                       fill = replicate_colors,
                                       filename = NULL)
  cosmic_variant_venns[[samp]] <- venn_plot_cosmic_obj
  
  variant_acast_cosmic[!is.na(variant_acast_cosmic)] <- 1
  variant_acast_cosmic <- data.frame(variant_acast_cosmic)
  variant_cosmic_common <- variant_acast_cosmic %>%
    dplyr::mutate(replicate_exist = paste0(L001, L002, L003, L004))
  
  variant_cosmic_common$id <- row.names(variant_acast_cosmic)
  variant_cosmic_common <-
    variant_cosmic_common[variant_cosmic_common$replicate_exist == "1111", ]
  
  common_var <- sample_variant_df$variant_id %in% variant_cosmic_common$id
  sample_subset_df <- sample_variant_df[common_var, ]
  sample_subset_df <- sample_subset_df[sample_subset_df$replicate == "L001", ]
  all_var_df <- rbind(all_var_df, sample_subset_df)
}

# Build a grid arrange gTree to visualize sequential Venn Diagrams
build_all <- ""
build_cosmic <- ""
for (venn_idx in 1:length(all_variant_venns)) {
  if (venn_idx == 1) {
    build_all <- "grid.arrange(gTree(children = all_variant_venns[[1]]),"
    build_cosmic <- "grid.arrange(gTree(children = cosmic_variant_venns[[1]]),"
  } else if (venn_idx < length(all_variant_venns)) {
    venn_ <- paste0("gTree(children = all_variant_venns[[", venn_idx, "]]),")
    venn_c <- paste0("gTree(children = cosmic_variant_venns[[", venn_idx,
                     "]]),")
    build_all <- paste(build_all, venn_)
    build_cosmic <- paste(build_cosmic, venn_c)
  } else {
    venn_ <- paste0("gTree(children = all_variant_venns[[", venn_idx,
                    "]]), ncol=6)")
    venn_c <- paste0("gTree(children = cosmic_variant_venns[[", venn_idx,
                     "]]), ncol=6)")
    build_all <- paste(build_all, venn_)
    build_cosmic <- paste(build_cosmic, venn_c)
  }
}

all_arrange_file <- file.path("figures", "venn_all_annotated_variants.pdf")
all_arrange_plot <- eval(parse(text = build_all))
ggplot2::ggsave(filename = all_arrange_file,
                plot = all_arrange_plot, height = 13, width = 13)

cosmic_file <- file.path("figures", "venn_cosmic_annotated_variants.pdf")
cosmic_arrange_plot <- eval(parse(text = build_cosmic))
ggplot2::ggsave(filename = cosmic_file,
                plot = cosmic_arrange_plot, height = 13, width = 13)

# Add groups together manually
all_var_df[all_var_df$sample_id == "KS25", "sample_id"] <- "KS25_primary"
all_var_df[all_var_df$sample_id == "KS26", "sample_id"] <- "KS26_F0"
all_var_df[all_var_df$sample_id == "KS27", "sample_id"] <- "KS27_F0"
all_var_df[all_var_df$sample_id == "KS28", "sample_id"] <- "KS28_F3"

all_var_df[all_var_df$base_id == "KS25", "base_id"] <- "KS25_KS26"
all_var_df[all_var_df$base_id == "KS26", "base_id"] <- "KS25_KS26"
all_var_df[all_var_df$base_id == "KS27", "base_id"] <- "KS27_KS28"
all_var_df[all_var_df$base_id == "KS28", "base_id"] <- "KS27_KS28"

# Write out the results
var_file <- file.path("results", "all_common_replicate_COSMIC_variants.tsv")
write.table(all_var_df, file = var_file, sep = "\t", row.names = FALSE)

# Assess consistency across sample passage
group_venns <- list()
for (sample_group in unique(all_var_df$base_id)) {
  all_called_var_subset <- all_var_df[all_var_df$base_id == sample_group, ]
  variant_acast_full_group <- reshape2::acast(all_called_var_subset,
                                              variant_id ~ sample_id,
                                              value.var = "variant_id")
  variant_list_full_group <- lapply(as.data.frame(variant_acast_full_group), 
                                    function(x) x[!is.na(x)])
  
  venn_colors <- replicate_colors[1:length(variant_list_full_group)]

  venn_plot_full_obj <- venn.diagram(x = variant_list_full_group,
                                     fill = venn_colors,
                                     margin = 0.16,
                                     cat.dist = 0.15,
                                     filename = NULL)
  group_venns[[sample_group]] <- venn_plot_full_obj
}

# Output visualization of group venn diagrams across passages
build_arrange_group <- ""
for (venn_idx in 1:length(group_venns)) {
  if (venn_idx == 1) {
    build_arrange_group <- "grid.arrange(gTree(children = group_venns[[1]]),"
  } else if (venn_idx < length(group_venns)) {
    venn_call <- paste0("gTree(children = group_venns[[", venn_idx, "]]),")
    build_arrange_group <- paste(build_arrange_group, venn_call)
  } else {
    venn_call <- paste0("gTree(children = group_venns[[", venn_idx,
                        "]]), ncol=3)")
    build_arrange_group <- paste(build_arrange_group, venn_call)
  }
}

group_file <- file.path("figures", "venn_group_annotated_variants.pdf")
group_arrange_plot <- eval(parse(text = build_arrange_group))
ggplot2::ggsave(filename = group_file,
                plot = group_arrange_plot, height = 10, width = 8)
