setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(Seurat)

suppressMessages({
  # extrafont::font_import()
  extrafont::loadfonts(device="postscript")
})

set.seed(42)
cells.sce = readRDS("../data/processed/cells.sce.rds")

# FigS1a
as.data.frame(colData(cells.sce)[,c("nCount_RNA","nFeature_RNA","discard")]) %>%
  ggplot(mapping = aes(x=nCount_RNA, y=nFeature_RNA)) +
  geom_point(aes(color = discard), size = .8) +
  scale_color_manual(values = c("#1F78B4","#E31A1C"), labels = c("Kept","Removed"), name = "") +
  scale_x_continuous(expand = c(0.15,0.1)) +
  xlab("UMI counts") +
  ylab("Number of Genes Recovered") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18),
    legend.title = element_blank(),
    legend.box.background = element_blank(),
    legend.position = c(.95, .05),
    legend.justification = c("right", "bottom"),
    legend.text = element_text(size = 15),
    legend.margin = margin(1, 5, 5, 5)
  )
ggsave(filename="sce_scatter.eps", plot=last_plot(), device="eps", path="../figures/FigS1/", width=5, height=4, dpi=300, family="Arial")

# FigS1b
as.data.frame(colData(cells.sce)[,c("subsets_Mito_percent","stim","high_subsets_Mito_percent")]) |>
  dplyr::mutate(stim = factor(stim, 
                              levels = c("Ctrl","LPS 500ug 3h","LPS 500ug 6h","LPS 500ug 1d","LPS 500ug 2d","LPS 10mg 3h",
                                         "LPS 10mg 6h","LPS 50mg 3h","LPS 50mg 6h","LPS >3w",
                                         "Poly(i:c) 10mg 3h","Poly(i:c) 10mg 6h","Poly(i:c) 20mg 3h","Poly(i:c) 20mg 6h","Poly(i:c) >3w","TNFalpha 500ug 6h"))) |>
  # dplyr::mutate(stim = dplyr::recode(stim, `LPS 500ug 3h`="LPS 500μg 3h",`LPS 500ug 6h`="LPS 500μg 6h",`LPS 500ug 1d`="LPS 500μg 1d",`LPS 500ug 2d`="LPS 500μg 2d",
  #                                    `TNFalpha 500ug 6h`="TNFalpha 500μg 6h")) |>  # grid converions failure
  ggplot(mapping = aes(x=stim, y=subsets_Mito_percent)) +
  geom_violin() +
  geom_jitter(mapping = aes(color = high_subsets_Mito_percent), size = 0.3, height = 0, width = .1) +
  scale_color_manual(values = c("#1F78B4","#E31A1C")) +
  theme_linedraw() +
  xlab("Experimental Conditions") +
  ylab("Percent Mitochondrial (%)") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)
  )
ggsave(filename="mito_violin.eps", plot=last_plot(), device="eps", path="../figures/FigS1/", width=6, height=5, dpi=300, family="Arial")

# FigS1c
as.data.frame(colData(cells.sce)[,c("subsets_Ribo_percent","stim","high_subsets_Ribo_percent")]) |>
  dplyr::mutate(stim = factor(stim, 
                              levels = c("Ctrl","LPS 500ug 3h","LPS 500ug 6h","LPS 500ug 1d","LPS 500ug 2d","LPS 10mg 3h",
                                         "LPS 10mg 6h","LPS 50mg 3h","LPS 50mg 6h","LPS >3w",
                                         "Poly(i:c) 10mg 3h","Poly(i:c) 10mg 6h","Poly(i:c) 20mg 3h","Poly(i:c) 20mg 6h","Poly(i:c) >3w","TNFalpha 500ug 6h"))) |>
  ggplot(mapping = aes(x=stim, y=subsets_Ribo_percent)) +
  geom_violin() +
  geom_jitter(mapping = aes(color = high_subsets_Ribo_percent), size = 0.3, height = 0, width = .1) +
  scale_color_manual(values = c("#1F78B4","#E31A1C")) +
  xlab("Experimental Conditions") +
  ylab("Percent Ribosomal (%)") +
  theme_linedraw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)
  )
ggsave(filename="ribo_violin.eps", plot=last_plot(), device="eps", path="../figures/FigS1/", width=6, height=5, dpi=300, family="Arial")


cells_postprocessed <- SeuratDisk::LoadH5Seurat("../data/processed/cells_postprocessed.h5Seurat", verbose=F)
# FigS1d
cells_postprocessed |> 
  FetchData(vars = c("cell_type_refined","treat")) |>
  dplyr::mutate(cell_type_refined = factor(cell_type_refined,
                                           levels = c("Somatotropes", "Lactotropes", "Corticotropes", "Melanotropes",
                                                      "Gonadotropes", "Thyrotropes", "Pou1f1 Progenitors", "Stem Cells",
                                                      "White Blood Cells", "Red Blood Cells", "Endothelial Cells",
                                                      "Pericytes", "Pituicytes", "Ambiguous"))) |>
  dplyr::mutate(treat = factor(treat,
                               levels = c("Saline","LPS","Poly(i:c)","TNFalpha"))) |>
  dplyr::group_by(treat, cell_type_refined) |>
  dplyr::summarise(n = n()) |>
  mutate(prop = n/sum(n)*100) |>
  ggplot(mapping = aes(fill=cell_type_refined, y=prop, x=treat)) +
  geom_bar(position="fill", stat = "identity", color = "white", size = .2, width = .7) +
  ylab("Proportion (%)") +
  scale_fill_manual(name = "Cell types", values = scales::hue_pal()(14)) +
  # scale_x_discrete(labels = c("WT\nPBS","WT\nLPS","GDKO\nPBS","GDKO\nLPS")) +
  scale_y_continuous(labels = scales::percent) +
  theme_linedraw() +
  theme(
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 16, angle = 90, hjust = 1, vjust = .5),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16)
  )
ggsave(filename="cell_proportion.eps", plot=last_plot(), device="eps", path="../figures/FigS1/", width=6, height=5, dpi=300, family="Arial")
