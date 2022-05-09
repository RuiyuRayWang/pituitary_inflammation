library(Seurat)
library(SeuratDisk)
library(Augur)

hpcs.lps <- LoadH5Seurat(file = "../data/processed/hpcs_lps_state_marked.h5Seurat", verbose = F)

hpcs.augur <- calculate_auc(hpcs.lps, cell_type_col = "cell_type_brief", label_col = "state", n_threads = 32)

saveRDS(hpcs.augur, file = "../data/processed/hpcs.augur.rds")

# https://gist.github.com/ramhiser/93fe37be439c480dc26c4bed8aab03dd
# https://www.r-graph-gallery.com/267-reorder-a-variable-in-ggplot2.html
augur.df <- hpcs.augur$AUC
df.tmp <- data.frame(cell_type = c("Som","Lac","Cort","Mel","Gonad","Thyro"), 
                     n_degs = c(sum(de.lps.markers.all.df$cell_type == "Som"),
                                sum(de.lps.markers.all.df$cell_type == "Lac"),
                                sum(de.lps.markers.all.df$cell_type == "Cort"),
                                sum(de.lps.markers.all.df$cell_type == "Mel"),
                                sum(de.lps.markers.all.df$cell_type == "Gonad"),
                                sum(de.lps.markers.all.df$cell_type == "Thyro")))
augur.df <- dplyr::left_join(augur.df, df.tmp, by = "cell_type")
augur.df <- augur.df %>% mutate_if(sapply(augur.df, is.character), as.factor)

my_color_palette <- scales::hue_pal()(length(unique(augur.df$cell_type)))
scales::show_col(my_color_palette)
col_pal.df <- as.data.frame(
  cbind(
    my_color_palette,
    c("Som","Cort","Lac","Gonad","Thyro")
  ))
colnames(col_pal.df) <- c("col_pal","cell_type")
augur.df <- dplyr::right_join(augur.df, col_pal.df, by = "cell_type")
augur.df <- dplyr::mutate(augur.df, cell_type = fct_reorder(cell_type, desc(auc)))