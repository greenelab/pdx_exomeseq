
library(UpSetR)
library(reshape2)
library(dplyr)

cosmic_file <- file.path('results', 'all_cosmic_variants.tsv')
cosmic_df <- readr::read_tsv(cosmic_file)
head(cosmic_df)

a_samples <- c('KS27', 'KS28', 'KS29', 'KS30')
patient_a <- cosmic_df %>% dplyr::filter(sample_name %in% a_samples)

a_melted <- reshape2::melt(patient_a, id.vars = 'sample_name', measure.vars = 'Gene.refGene')
a_pivot <- reshape2::dcast(a_melted, value ~ sample_name, fun.aggregate = function(x){ length(x) })
a_pivot[a_pivot == 2] <- 1

upset(a_pivot, order.by = 'freq', sets = a_samples, keep.order = TRUE,
      queries = list(list(query = intersects, params = list(a_samples),
                          color='orange', active=T)), mb.ratio = c(0.7, 0.3),
      sets.bar.color = c("#5eff2d", "#5eff2d", "#e8257c", "#e8257c"),
      text.scale = c(2, 1.5, 2, 1.5, 2, 1.2))

b_samples <- c('008-F0', '008-F5', '018-F0', '018-F5', '019-F0', '019-F5', 'KS25', 'KS26')
patient_b <- cosmic_df %>% dplyr::filter(sample_name %in% b_samples)

b_melted <- reshape2::melt(patient_b, id.vars = 'sample_name', measure.vars = 'Gene.refGene')
b_pivot <- reshape2::dcast(b_melted, value ~ sample_name, fun.aggregate = function(x){ length(x) })
b_pivot[b_pivot == 2] <- 1

upset(b_pivot, order.by = 'freq', sets = b_samples, keep.order = TRUE,
      queries = list(list(query = intersects, params = list(b_samples),
                          color='orange', active=T)), mb.ratio = c(0.7, 0.3),
      sets.bar.color = c("#5eff2d", "#5eff2d", "#332dff", "#332dff",
                         "#e8257c", "#e8257c", "#32d2ff", "#32d2ff"),
      text.scale = c(2, 1.5, 2, 1.5, 2, 1.2))

upset_a_file <- file.path('figures', 'upset_sample_a.pdf')
pdf(upset_a_file, height = 6, width = 7)
upset(a_pivot, order.by = 'freq', sets = a_samples, keep.order = TRUE,
      queries = list(list(query = intersects, params = list(a_samples),
                          color='orange', active=T)), mb.ratio = c(0.7, 0.3),
      sets.bar.color = c("#5eff2d", "#5eff2d", "#e8257c", "#e8257c"),
      text.scale = c(2, 1.5, 2, 1.5, 2, 1.2))
dev.off()
dev.off()

upset_b_file <- file.path('figures', 'upset_sample_b.pdf')
pdf(upset_b_file, height = 6, width = 7)
upset(b_pivot, order.by = 'freq', sets = b_samples, keep.order = TRUE,
      queries = list(list(query = intersects, params = list(b_samples),
                          color='orange', active=T)), mb.ratio = c(0.7, 0.3),
      sets.bar.color = c("#5eff2d", "#5eff2d", "#332dff", "#332dff",
                         "#e8257c", "#e8257c", "#32d2ff", "#32d2ff"),
      text.scale = c(2, 1.5, 2, 1.5, 2, 1.2))
dev.off()
dev.off()
