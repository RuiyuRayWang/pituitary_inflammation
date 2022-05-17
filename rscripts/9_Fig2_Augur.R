library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(ggrepel)
library(Augur)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

suppressMessages(
  extrafont::loadfonts(device="postscript")
)

hpcs.lps <- LoadH5Seurat(file = "../data/processed/hpcs_lps_state_marked.h5Seurat", verbose = F)

test_method = "MAST"
de.state.mks <- read.csv(file = paste0("../outs/de_state_markers_",test_method,".csv"))
# de.state.mks.uniq <- de.state.mks[!duplicated(de.state.mks$gene),]

hpcs.augur <- calculate_auc(hpcs.lps, cell_type_col = "cell_type_brief", label_col = "state", n_threads = 32)

saveRDS(hpcs.augur, file = "../data/processed/hpcs.augur.rds")

# https://gist.github.com/ramhiser/93fe37be439c480dc26c4bed8aab03dd
# https://www.r-graph-gallery.com/267-reorder-a-variable-in-ggplot2.html
hpcs.augur <- readRDS("../data/processed/hpcs.augur.rds")
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
ggplot(augur.df, mapping = aes(x = cell_type, y = auc)) +
  coord_flip() +
  geom_point(aes(fill=cell_type), size=6, shape=21, stroke=0) +
  geom_segment(aes(x=cell_type, xend=cell_type, y=0, yend=auc, color = cell_type), size = 2) +
  scale_fill_manual(values = augur.df$col_pal,
                    breaks = augur.df$cell_type) +
  scale_color_manual(values = augur.df$col_pal,
                     breaks = augur.df$cell_type) +
  geom_text(mapping = aes(label = round(auc, digits = 2)), 
            vjust = -.75,
            # hjust = -0.25,
            size = 6) +
  # ylim(c(-0.02, 1.2)) +
  scale_y_continuous(limits = c(0.0, 1.05),expand = c(0, 0)) +  # Force Y-axis start from 0
  xlab(label = "Cell type") +
  ylab(label = "Augur AUC") +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(size = 2, fill = NA),
    # axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 20),
    legend.position = "none"
  )
ggsave(
  filename = "augur_lolipop.eps",
  device = "eps", 
  plot = last_plot(), 
  path = "../figures/Fig2/", 
  dpi = 300, width = 4, height = 5,
  family = "Arial"
)

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
    panel.border = element_rect(size = 1.6, fill = NA),
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

