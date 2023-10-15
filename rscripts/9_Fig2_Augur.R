setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(Augur)
library(Seurat)
library(tidyverse)
library(viridis)
source("plot_umap_refactored.R")

suppressMessages({
  # extrafont::font_import()
  extrafont::loadfonts(device="postscript")
})

hpcs.augur <- readRDS(file = "../data/processed/hpcs.augur.rds")

# ======================================================================================================
# Augur AUC on UMAP
library(SeuratDisk)
hpcs.lps <- LoadH5Seurat("../data/processed/hpcs_lps_state_marked.h5Seurat")
hpcs.lps$cell_type = hpcs.lps$cell_type_brief
hpcs.lps$barcode <- NULL

u <- plot_umap_refactored(
  augur = hpcs.augur, sc = hpcs.lps, palette = "OrRd", reduction = "umap", cell_type_col = "cell_type", 
  size_sm = 14, size_lg = 16, lgd.pos = c(.12,.05), lgd.name = "Augur\nAUC", direction = 1
)
ggsave(filename = "hpcs_augurscore_umap.eps", plot = u, device = "eps", path = "../figures/Fig2/", width = 4, height = 4, dpi = 300, family = "Arial")

# ======================================================================================================
# Lollipop plot
aucs = hpcs.augur$AUC
size_sm = 12
size_lg = 14
range = range(aucs$auc)
expand = abs(diff(range)) * 0.25
p = aucs |> ggplot(aes(x = reorder(cell_type, auc), y = auc)) + 
  geom_hline(aes(yintercept = 0.5), linetype = "dotted", size = 0.5) + 
  geom_point(aes(color = cell_type), size = 2) + 
  geom_text(aes(label = format(auc, digits = 2), y = ifelse(auc < 0.5, 0.5, auc)), size = 5, nudge_y = expand, hjust = 0.5, color = "black") + 
  geom_segment(aes(color = cell_type, xend = cell_type, yend = 0.5)) + scale_y_continuous("Augur AUC", limits = c(min(range[1] - expand, 0.5), range[2] + expand * 1.5)) + 
  scale_color_manual(values = scales::hue_pal()(13)[c(9,10,2,6,1,11)]) +
  coord_flip() + 
  theme_bw() + 
  theme(axis.title.x = element_text(size = size_lg), axis.text.x = element_text(size = size_sm, color = "black"), 
        axis.title.y = element_blank(), axis.text.y = element_text(size = size_lg-0.6, color = "black"),
        panel.grid = element_blank(), 
        strip.text = element_text(size = size_lg), strip.background = element_blank(), 
        axis.line.y = element_blank(), axis.line.x = element_blank(), 
        legend.position = "none", legend.text = element_text(size = size_sm), 
        legend.title = element_blank(), legend.key.size = unit(0.6, "lines"), 
        legend.margin = margin(rep(0, 4)), legend.background = element_blank(), 
        plot.title = element_text(size = size_lg, hjust = 0.5))
ggsave(
  filename = "augur_lolipop.eps",
  device = "eps", 
  plot = last_plot(), 
  path = "../figures/Fig2/", 
  dpi = 300, width = 3, height = 3.6,
  family = "Arial"
)

# ======================================================================================================

test_method = "MAST"
de.state.mks <- read.csv(file = paste0("../outs/de_state_markers_",test_method,".csv"))

augur.df <- hpcs.augur$AUC
df.tmp <- data.frame(cell_type = c("Som","Lac","Cort","Mel","Gonad","Thyro"), 
                     n_degs = c(sum(de.state.mks$cell_type == "Som"),
                                sum(de.state.mks$cell_type == "Lac"),
                                sum(de.state.mks$cell_type == "Cort"),
                                sum(de.state.mks$cell_type == "Mel"),
                                sum(de.state.mks$cell_type == "Gonad"),
                                sum(de.state.mks$cell_type == "Thyro")))
augur.df <- dplyr::left_join(augur.df, df.tmp, by = "cell_type")
augur.df <- augur.df %>% mutate_if(sapply(augur.df, is.character), as.factor)

my_color_palette <- scales::hue_pal()(13)[c(1,2,9,6,10,13)]
scales::show_col(my_color_palette)
col_pal.df <- as.data.frame(
  cbind(
    my_color_palette,
    c("Som","Lac","Cort","Mel","Gonad","Thyro")
  ))
colnames(col_pal.df) <- c("col_pal","cell_type")
augur.df <- dplyr::right_join(augur.df, col_pal.df, by = "cell_type")
augur.df <- dplyr::mutate(augur.df, cell_type = fct_reorder(cell_type, dplyr::desc(auc)))

augur.df$cell_type <- factor(augur.df$cell_type, levels = rev(c("Cort","Som","Mel","Lac","Gonad","Thyro")))

pearson.res <- cor.test(augur.df$n_degs, augur.df$auc, method = "pearson", conf.level = 0.95)
ggplot(augur.df, aes(x = n_degs, y = auc)) +
  geom_point(aes(color = cell_type), size=5) +
  geom_smooth(method="lm", se = FALSE, linetype = "dashed") +
  # geom_text_repel(aes(label = cell_type), box.padding = 1, size = 4, force = 5) +
  geom_label_repel(data = subset(augur.df, cell_type %in% c("Som")), aes(label = cell_type), 
                   nudge_x = -30, nudge_y = 0.04, box.padding = 0.25, point.padding = 0.4, size = 5, force = 5) +
  geom_label_repel(data = subset(augur.df, cell_type %in% c("Lac")), aes(label = cell_type),
                   nudge_x = 100,nudge_y = 0.005, box.padding = 0.25, point.padding = 0.4, size = 5, force = 5) +
  geom_label_repel(data = subset(augur.df, cell_type %in% c("Cort")), aes(label = cell_type), 
                   nudge_x = -30, nudge_y = -0.05, box.padding = 0.25, point.padding = 0.4, size = 5, force = 5) +
  geom_label_repel(data = subset(augur.df, cell_type %in% c("Mel")), aes(label = cell_type),
                   nudge_x = 100,nudge_y = -0.01, box.padding = 0.25, point.padding = 0.4, size = 5, force = 5) +
  geom_label_repel(data = subset(augur.df, cell_type %in% c("Gonad")), aes(label = cell_type), 
                   nudge_x = 120, nudge_y = 0, box.padding = 0.25, point.padding = 0.4, size = 5, force = 1) +
  geom_label_repel(data = subset(augur.df, cell_type %in% c("Thyro")), aes(label = cell_type), 
                   nudge_y = 0.04, box.padding = 0.25, point.padding = 0.4, size = 5, force = 5) +
  xlab("Number of DE genes") +
  ylab("Augur AUC") +
  scale_color_manual(values = augur.df$col_pal,
                     breaks = augur.df$cell_type) +
  scale_x_continuous(minor_breaks = seq(0,500,100)) +
  scale_y_continuous(minor_breaks = seq(0,0.9,0.05)) +
  annotate("text", x = 300, y = 0.92, size = 6,
           label = paste0("R = ",round(pearson.res$estimate, digits = 3),", p = ", sprintf("%.5f",pearson.res$p.value)), 
           color = "blue") +
  theme(
    panel.border = element_rect(size = .6, fill = NA),
    panel.background = element_blank(),
    panel.grid = element_line(color = "#CCCCCC", size = .6),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.position = "none"
  )

ggsave(
  filename = "augur_auc_de_scatter.eps",
  device = "eps", 
  plot = last_plot(), 
  path = "../figures/Fig2/", 
  dpi = 300, width = 4.5, height = 4.5,
  family = "Arial"
)
