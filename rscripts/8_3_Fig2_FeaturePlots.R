library(Seurat)
library(SeuratDisk)
library(RColorBrewer)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

suppressMessages(
  extrafont::loadfonts(device="postscript")
)

hpcs.lps <- LoadH5Seurat("../data/processed/hpcs_lps_state_marked.h5Seurat")

DimPlot(hpcs.lps, reduction = "umap", group.by = "cell_type_brief", label = T) +
  NoAxes() +
  ggtitle("") +
  theme(
    text = element_text(size = 16)
  )
ggsave(
  filename = "umap_cell_type.eps",
  plot = last_plot(),
  device = "eps", 
  path = "../figures/Fig2/umap_featureplots_yt/",
  width = 5, height = 5, 
  dpi = 300
)

go0051047 <- c(
  "Rab27b","Nr3c1","Rtn4","S100a8","Cartpt","S100a9","Apbb1","Crhr1","Ffar4","Sphk1","Hif1a","Baiap3","Tac1","Irs2","Dpysl2","Pcsk1"  ## GO:0051047 positive regulation of secretion
  )
for (feat in go0051047){
  FeaturePlot(hpcs.lps, features = feat, reduction = "umap") +
    scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd")) +
    ggtitle("") +
    NoAxes() +
    annotate(label = feat, geom = "text", x = -11, y = -12, hjust = 1, size = 6) +
    theme(
      legend.position = c(0.86, 0.16),
      panel.border = element_rect(colour = "black", fill=NA, size=.5)
    )
  ggsave(
    filename = paste0(feat,"_umap.eps"),
    plot = last_plot(),
    device = "eps",
    path = "../figures/Fig2/umap_featureplots_yt/go0051047",
    width = 5, height = 5,
    dpi = 300,
    family = "Arial"
  )
}

go0033673 <- c(
  "Dusp26","Trib1","Gadd45a","Rgs2","Dusp1","Gprc5a","Pik3ip1","Fabp4","Irs2","Gstp1","Ptpn2"  ## GO:0033673 negative regulation of kinase activity
)
for (feat in go0033673){
  FeaturePlot(hpcs.lps, features = feat, reduction = "umap") +
    scale_colour_gradientn(colours = brewer.pal(n = 9, name = "YlOrRd")) +
    ggtitle("") +
    NoAxes() +
    annotate(label = feat, geom = "text", x = -11, y = -12, hjust = 1, size = 6) +
    theme(
      legend.position = c(0.86, 0.16),
      panel.border = element_rect(colour = "black", fill=NA, size=.5)
    )
  ggsave(
    filename = paste0(feat,"_umap.eps"),
    plot = last_plot(),
    device = "eps",
    path = "../figures/Fig2/umap_featureplots_yt/go0033673",
    width = 5, height = 5,
    dpi = 300,
    family = "Arial"
  )
}