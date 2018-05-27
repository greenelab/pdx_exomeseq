
library(UpSetR)
library(reshape2)
library(dplyr)

cosmic_file <- file.path('results', 'all_cosmic_variants.tsv')
cosmic_df <- readr::read_tsv(cosmic_file) %>%
    dplyr::mutate(base_sample_id = sapply(final_id, function(x) unlist(strsplit(x, '-'))[1]))
head(cosmic_df)

for (sample_group in unique(cosmic_df$base_sample_id)) {
    
    upset_fig_file <- file.path('figures', 'upset', paste0('upset_sample_', sample_group, '.pdf'))
   
    
    patient_df <- cosmic_df %>% dplyr::filter(base_sample_id == sample_group)
    
    sample_set <- sort(unique(patient_df$final_id), decreasing = TRUE)

    patient_df_melt <- reshape2::melt(patient_df, id.vars = 'final_id', measure.vars = 'Gene.refGene')
    patient_pivot <- reshape2::dcast(patient_df_melt, value ~ final_id, fun.aggregate = function(x) length(x) )
    patient_pivot[patient_pivot == 2] <- 1
    
    pdf(upset_fig_file, height = 6, width = 7, onefile = FALSE)         
    upset(patient_pivot, order.by = 'freq', sets = sample_set, keep.order = TRUE,
          queries = list(list(query = intersects, params = sample_set,
                              color = 'orange', active = T)), mb.ratio = c(0.7, 0.3),
          text.scale = c(1.8, 1.5, 1.8, 1.5, 1.8, 1.2))
    dev.off()
}

file = file.path('results', 'all_cosmic_prefiltered_variants.tsv')
prefiltered_cosmic_df <- readr::read_tsv(file) %>%
    dplyr::mutate(base_sample_id = sapply(final_id, function(x) unlist(strsplit(x, '-'))[1]))
head(prefiltered_cosmic_df)

for (sample_group in unique(prefiltered_cosmic_df$base_sample_id)) {
    
    upset_fig_file <- file.path('figures', 'upset', 'prefiltered',
                                paste0('upset_sample_', sample_group, '.pdf'))
   
    
    patient_df <- prefiltered_cosmic_df %>% dplyr::filter(base_sample_id == sample_group)
    
    sample_set <- sort(unique(patient_df$final_id), decreasing = TRUE)

    patient_df_melt <- reshape2::melt(patient_df, id.vars = 'final_id', measure.vars = 'Gene.refGene')
    patient_pivot <- reshape2::dcast(patient_df_melt, value ~ final_id, fun.aggregate = function(x) length(x) )
    patient_pivot[patient_pivot >= 2] <- 1
    
    pdf(upset_fig_file, height = 6, width = 7, onefile = FALSE)         
    upset(patient_pivot, order.by = 'freq', sets = sample_set, keep.order = TRUE,
          queries = list(list(query = intersects, params = sample_set,
                              color = 'orange', active = T)), mb.ratio = c(0.7, 0.3),
          text.scale = c(1.8, 1.5, 1.8, 1.5, 1.8, 1.2))
    dev.off()
}
