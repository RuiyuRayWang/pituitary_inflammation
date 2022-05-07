library(Seurat)
library(SeuratDisk)
library(scales)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

suppressMessages(
  extrafont::loadfonts(device="postscript")
)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

cells <- LoadH5Seurat("../data/cells_postprocessed.h5Seurat", verbose = F)

hpcs.lps <- LoadH5Seurat("../data/hpcs_lps_state_marked.h5Seurat", verbose = F)

# Fig1b
## Reorder factor levels
cells$cell_type_brief <- factor(cells$cell_type_brief, levels = c("Som", "Lac", "Cort", "Mel",
                                                                  "Gonad", "Thyro", "Pou1f1", "Stem",
                                                                  "WBCs", "RBCs", "Endo",
                                                                  "Peri", "Pitui", "Ambig"))

DimPlot(cells, reduction = "umap.int", cols = c(hue_pal()(13)[c(1,2,9,6,10,13,3,4,5,7,8,9,11)],"#CCCCCC")) + 
  NoLegend() +
  xlab("UMAP_1") + ylab("UMAP_2") +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 18)
  )
ggsave(
  filename = "umapint_celltypebrief_clean.eps", 
  plot = last_plot(), 
  device = "eps", 
  path = "../figures/Fig1", 
  dpi = 300,
  family = "Arial"
  )

# Fig1d
c_features <- c("Gh", "Ghrhr", "Pomc", "Pax7", "Prl", "Cga", "Lhb", "Fshb", "Tshb", "Trhr")

for (feat in c_features){
  FeaturePlot(cells, features = feat, reduction = "umap.int") +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu"))) +
    ggtitle("") +
    NoAxes() +
    annotate(label = feat, geom = "text", x = 13, y = 10, hjust = 1, size = 7) +
    theme(
      legend.position = c(0.06, 0.24),
      panel.border = element_rect(colour = "black", fill=NA, size=.5)
    )
  ggsave(
    filename = paste0(feat,"_umapint.eps"),
    plot = last_plot(),
    device = "eps",
    path = "../figures/Fig1/F1d/", 
    width = 4, height = 4,
    dpi = 300,
    family = "Arial"
  )
}

## DEPRECATED
# f.list <- FeaturePlot(cells, features = c_features, reduction = "umap.int", combine = FALSE)
# for (i in 1:length(f.list)){
#   f.list[[i]] <- f.list[[i]] + 
#     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu"))) +
#     ggtitle("") + 
#     NoAxes() + 
#     annotate(label = c_features[i], geom = "text", x = 13, y = 10.8, hjust = 1, size = 6) +
#     theme(
#       legend.position = c(0.06, 0.24),
#       panel.border = element_rect(colour = "black", fill=NA, size=.5)
#     )
# }
# wrap_plots(f.list) + plot_layout(ncol = 5, nrow = 2)

# Fig1e
DefaultAssay(hpcs.lps) <- "RNA"
n_epochs = 400
a = 0.88; b = 1.1
hpcs.lps <- RunUMAP(hpcs.lps, reduction = "pca", dims = 1:40, n.epochs = n_epochs, a = a, b = b)

DimPlot(hpcs.lps, group.by = "state", reduction = "umap", cols = brewer.pal(n = 12, name = "Paired")[c(8,2)], pt.size = .6) +
  NoAxes() +
  NoLegend() +
  ggtitle(NULL)
# FetchData(
#   object = hpcs.lps, 
#   vars = c("UMAP_1", "UMAP_2", "state", "cell_type_brief")
# ) %>%
#   mutate(state = factor(state, levels = c("Healthy","Inflammation")),
#          cell_type_brief = factor(cell_type_brief, levels = c("Som","Lac","Cort","Mel","Gonad","Thyro"))) %>%
#   ggplot(aes(x = UMAP_1, y = UMAP_2)) +
#   geom_point(aes(color = cell_type_brief, shape = state, fill = state), size = 1.2) +
#   scale_color_manual(values = hue_pal()(13)[c(1,2,9,6,10,11)]) +
#   scale_shape_manual(values = c(21,23)) +
#   scale_fill_manual(values = brewer.pal(n = 8, name = "Paired")[c(1,8)]) +
#   theme(
#     axis.line = element_blank(),
#     axis.title = element_blank(),
#     axis.text = element_blank(),
#     axis.ticks = element_blank(),
#     panel.background = element_blank(),
#     legend.position = "none",
#   )
ggsave(
  filename = "umap_hpcs_lps_clean.eps",
  plot = last_plot(), 
  device = "eps", 
  path = "../figures/Fig1/", 
  width = 4, height = 4,
  dpi = 300,
  )

# Fig1f
DefaultAssay(hpcs.lps) <- "RNA"
hpcs.lps$state <- factor(hpcs.lps$state, levels = c("Healthy", "Inflammation"))
## Som
som.lps <- subset(hpcs.lps, subset = cell_type_brief == "Som")
DefaultAssay(som.lps) <- "RNA"
som.lps <- som.lps %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:40)
DimPlot(som.lps, group.by = "state", reduction = "umap", cols = brewer.pal(n = 12, name = "Paired")[c(5,6)]) +
  NoAxes() +
  NoLegend() +
  ggtitle(NULL)
ggsave(
  filename = "umap_som_lps_clean.eps",
  plot = last_plot(), 
  device = "eps", 
  path = "../figures/Fig1/F1f/", 
  width = 5, height = 3,
  dpi = 300,
)
## Lac
lac.lps <- subset(hpcs.lps, subset = cell_type_brief == "Lac")
DefaultAssay(lac.lps) <- "RNA"
lac.lps <- lac.lps %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30)
DimPlot(lac.lps, group.by = "state", reduction = "umap", cols = brewer.pal(n = 12, name = "Paired")[c(7,8)]) +
  NoAxes() +
  NoLegend() +
  ggtitle(NULL)
ggsave(
  filename = "umap_lac_lps_clean.eps",
  plot = last_plot(), 
  device = "eps", 
  path = "../figures/Fig1/F1f/", 
  width = 4, height = 3,
  dpi = 300,
)
## Cort
cort.lps <- subset(hpcs.lps, subset = cell_type_brief == "Cort")
DefaultAssay(cort.lps) <- "RNA"
cort.lps <- cort.lps %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:25)
DimPlot(cort.lps, group.by = "state", reduction = "umap", cols = brewer.pal(n = 12, name = "Paired")[c(1,2)]) +
  NoAxes() +
  NoLegend() +
  ggtitle(NULL)
ggsave(
  filename = "umap_cort_lps_clean.eps",
  plot = last_plot(), 
  device = "eps", 
  path = "../figures/Fig1/F1f/", 
  width = 3, height = 3,
  dpi = 300,
)
## Mel
mel.lps <- subset(hpcs.lps, subset = cell_type_brief == "Mel")
DefaultAssay(mel.lps) <- "RNA"
mel.lps <- mel.lps %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30)
DimPlot(mel.lps, group.by = "state", reduction = "umap", cols = brewer.pal(n = 12, name = "Paired")[c(3,4)]) +
  NoAxes() +
  NoLegend() +
  ggtitle(NULL)
ggsave(
  filename = "umap_mel_lps_clean.eps",
  plot = last_plot(), 
  device = "eps", 
  path = "../figures/Fig1/F1f/", 
  width = 2, height = 3,
  dpi = 300,
)
## Gonad
gonad.lps <- subset(hpcs.lps, subset = cell_type_brief == "Gonad")
DefaultAssay(gonad.lps) <- "RNA"
gonad.lps <- gonad.lps %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30)
DimPlot(gonad.lps, group.by = "state", reduction = "umap", cols = c("#7570B3","#BEBADA")) +
  NoAxes() +
  NoLegend() +
  ggtitle(NULL)
ggsave(
  filename = "umap_gonad_lps_clean.eps",
  plot = last_plot(), 
  device = "eps", 
  path = "../figures/Fig1/F1f/", 
  width = 2, height = 3,
  dpi = 300,
)
## Thyro
thyro.lps <- subset(hpcs.lps, subset = cell_type_brief == "Thyro")
DefaultAssay(thyro.lps) <- "RNA"
thyro.lps <- thyro.lps %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:20)
DimPlot(thyro.lps, group.by = "state", reduction = "umap", cols = c("#F781BF","#E7298A")) +
  NoAxes() +
  NoLegend() +
  ggtitle(NULL)
ggsave(
  filename = "umap_thyro_lps_clean.eps",
  plot = last_plot(), 
  device = "eps", 
  path = "../figures/Fig1/F1f/", 
  width = 2, height = 2,
  dpi = 300,
)
