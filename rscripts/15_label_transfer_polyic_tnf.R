setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(Seurat)
# library(SeuratDisk)
library(tidyverse)
library(tidydr)
library(ggpubr)
library(RColorBrewer)

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
  scale_size_manual(values = c(1.2,1.2,.6,.6)) +
  theme_dr(xlength = .2, ylength = .2) +
  theme(
    legend.position = "none",
    axis.title = element_text(hjust = 0)
  )

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
  scale_size_manual(values = c(1, .3)) +
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

p3 <- DimPlot(
  query.hpcs, reduction = "ref.umap", group.by = "stim", label = TRUE,
  label.size = 3, label.box = T, repel = TRUE) + 
  ggtitle("Query treatments") +
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

ggarrange(p1,NULL,(p2/p3), ncol = 3, widths = c(1,.1,.5)) |>
  ggexport(filename = "../figures/polyic_tnf.pdf", width = 9, height = 5)
