setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(Seurat)
# library(SeuratDisk)
library(tidyverse)
library(tidydr)
library(ggpubr)
library(RColorBrewer)

suppressMessages(
  extrafont::loadfonts(device="postscript")
)

cells <- SeuratDisk::LoadH5Seurat("../data/processed/cells_postprocessed.h5Seurat")
hpcs.lps <- SeuratDisk::LoadH5Seurat("../data/processed/hpcs_lps_state_marked.h5Seurat", verbose = F)

query.hpcs <- subset(cells, subset = treat %in% c("Poly(i:c)","TNFalpha")) |>
  subset(subset = cell_type_brief %in% c("Som","Lac","Cort","Gonad","Mel","Thyro")) |>
  NormalizeData() |> FindVariableFeatures() |> ScaleData() |> RunPCA()

hpcs.lps$id = "Reference"
query.hpcs$id = "Query"

hpcs.anchors <- FindTransferAnchors(reference = hpcs.lps, query = query.hpcs, dims = 1:30, reference.reduction = "pca")
# predictions <- TransferData(anchorset = hpcs.anchors, refdata = hpcs.lps$state, dims = 1:30)
# query.hpcs <- AddMetaData(query.hpcs, metadata = predictions)

hpcs.lps <- RunUMAP(hpcs.lps, dims = 1:50, reduction = "pca", return.model = TRUE)
query.hpcs <- MapQuery(
  anchorset = hpcs.anchors, 
  reference = hpcs.lps, 
  query = query.hpcs,
  refdata = list(state = "state"), 
  reference.reduction = "pca", 
  reduction.model = "umap"
)
hpcs.lps$state <- factor(hpcs.lps$state, levels = c("Healthy","Inflammation"))
hpcs.lps$cell_type_brief <- factor(hpcs.lps$cell_type_brief, levels = c("Som", "Lac", "Cort", "Mel", "Gonad", "Thyro"))
query.hpcs$predicted.state <- factor(query.hpcs$predicted.state, levels = c("Healthy","Inflammation"))

p1 <- Embeddings(object = hpcs.lps, reduction = "umap") |>
  as.data.frame() |>
  rownames_to_column("cell_id") |>
  left_join(y = FetchData(object = hpcs.lps, vars = c("id","state")) |> rownames_to_column("cell_id"), by = "cell_id") |>
  column_to_rownames("cell_id") |>
  bind_rows(
    Embeddings(object = query.hpcs, reduction = "ref.umap") |>
      as.data.frame() |>
      rename(UMAP_1 = refUMAP_1, UMAP_2 = refUMAP_2) |>
      rownames_to_column("cell_id") |>
      left_join(y = FetchData(object = query.hpcs, vars = c("id","predicted.state")) |> rownames_to_column("cell_id"), by = "cell_id") |>
      rename(state = predicted.state) |>
      column_to_rownames("cell_id")
  ) |>
  unite(col = "id_state", id, state, sep = "_") |>
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = id_state, size = id_state)) +
  scale_color_manual(values = brewer.pal(12,"Paired")[c(2,8,1,7)]) +
  scale_size_manual(values = c(1,1,.5,.5)) +
  theme_dr(xlength = .16, ylength = .16) +
  theme(
    legend.position = "none",
    axis.title = element_text(hjust = 0)
  )
ggsave(filename = "p1.eps", plot = p1, device = "eps", path = "../figures/FigS2/", width = 6, height = 6, dpi = 300, family = "Arial")

hpcs.merged <- merge(hpcs.lps, query.hpcs)
p2 <- Embeddings(object = hpcs.lps, reduction = "umap") |>
  as.data.frame() |>
  rownames_to_column("cell_id") |>
  left_join(y = FetchData(object = hpcs.lps, vars = c("id","cell_type_brief")) |> rownames_to_column("cell_id"), by = "cell_id") |>
  column_to_rownames("cell_id") |>
  bind_rows(
    Embeddings(object = query.hpcs, reduction = "ref.umap") |>
      as.data.frame() |>
      rename(UMAP_1 = refUMAP_1, UMAP_2 = refUMAP_2) |>
      rownames_to_column("cell_id") |>
      left_join(y = FetchData(object = query.hpcs, vars = c("id","cell_type_brief")) |> rownames_to_column("cell_id"), by = "cell_id") |>
      column_to_rownames("cell_id")
  ) |>
  mutate(cell_type_brief = factor(cell_type_brief, levels = c("Som","Lac","Cort","Mel","Gonad","Thyro"))) |>
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = cell_type_brief, size = id)) +
  ggtitle("Cell types") +
  scale_color_manual(values = scales::hue_pal()(13)[c(1,2,9,6,10,13)]) +
  scale_size_manual(values = c(1, .2)) +
  theme_dr() +
  theme(
    legend.position = "none",
    panel.border = element_rect(fill = NA, color = "#696969"),
    axis.line.x.bottom = element_blank(),
    axis.line.y.left = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
p2 = LabelClusters(plot = p2, id = "cell_type_brief", box = T)
ggsave(filename = "p2.eps", plot = p2, device = "eps", path = "../figures/FigS2/", width = 3, height = 3, dpi = 300, family = "Arial")

p3 <- Embeddings(object = query.hpcs, reduction = "ref.umap") |>
  as.data.frame() |>
  rename(UMAP_1 = refUMAP_1, UMAP_2 = refUMAP_2) |>
  rownames_to_column("cell_id") |>
  left_join(y = FetchData(object = query.hpcs, vars = c("id","stim")) |> rownames_to_column("cell_id"), by = "cell_id") |>
  column_to_rownames("cell_id") |>
  mutate(labels = recode(stim, `Poly(i:c) >3w`="Poly(i:c) > 3w", `Poly(i:c) 10mg 3h`="Poly(i:c) 10mg/kg 3h", `Poly(i:c) 10mg 6h`="Poly(i:c) 10mg/kg 6h",
                         `Poly(i:c) 20mg 3h`="Poly(i:c) 20mg/kg 3h", `Poly(i:c) 20mg 6h`="Poly(i:c) 20mg/kg 6h", `TNFalpha 500ug 6h`="TNF-α 500ug/kg 6h")) |>
  mutate(labels = factor(labels, levels = c("Poly(i:c) > 3w","Poly(i:c) 10mg/kg 3h","Poly(i:c) 10mg/kg 6h","Poly(i:c) 20mg/kg 3h","Poly(i:c) 20mg/kg 6h","TNF-α 500ug/kg 6h"))) |>
  ggplot(mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = labels), size = .8) +
  ggtitle("Treatments") +
  scale_color_manual(values = brewer.pal(8,"Dark2")[1:6]) +
  theme_dr() +
  theme(
    legend.position = "none",
    panel.border = element_rect(fill = NA, color = "#696969"),
    axis.line.x.bottom = element_blank(),
    axis.line.y.left = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
p3 = LabelClusters(plot = p3, id = "labels", box = T, size = 2.5)
ggsave(filename = "p3.eps", plot = p3, device = "eps", path = "../figures/FigS2/", width = 3, height = 3, dpi = 300, family = "Arial")

ggarrange(p1,NULL,(p2/p3), ncol = 3, widths = c(1,.1,.5)) |>
  ggexport(filename = "../figures/FigS2/polyic_tnf.pdf", width = 9, height = 5)

p4 <- FetchData(object = query.hpcs, vars = c("stim","predicted.state")) |>
  dplyr::group_by(stim, predicted.state) |>
  dplyr::summarise(n = n()) |>
  mutate(labels = recode(stim, `Poly(i:c) >3w`="Poly(i:c) > 3w", `Poly(i:c) 10mg 3h`="Poly(i:c) 10mg/kg 3h", `Poly(i:c) 10mg 6h`="Poly(i:c) 10mg/kg 6h",
                         `Poly(i:c) 20mg 3h`="Poly(i:c) 20mg/kg 3h", `Poly(i:c) 20mg 6h`="Poly(i:c) 20mg/kg 6h", `TNFalpha 500ug 6h`="TNF-α 500ug/kg 6h")) |>
  ggplot(mapping = aes(fill=labels, y=n, x=predicted.state)) +
  geom_bar(stat = "identity", color = "white", size = .5, width = .5) +
  ylab("Numbre of cells") +
  xlab("Predicted state") +
  scale_fill_manual(name = "Treatment", values = brewer.pal(8,"Dark2")[1:6]) +
  theme_linedraw() +
  theme(
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    axis.title.x = element_text(size = 18),
    axis.text.x = element_text(size = 16),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(size = 16)
  )
ggsave(filename = "p4.eps", plot = p4, device = "eps", path = "../figures/FigS2/", width = 6, height = 5.4, dpi = 300, family = "Arial")

p4 |>
  ggexport(filename = "../figures/polyic_tnf_barproportion.pdf", width = 6, height = 6.4)