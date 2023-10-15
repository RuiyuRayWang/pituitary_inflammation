setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(scWidgets)
# library(clusterProfiler)
# library(org.Mm.eg.db)
# library(msigdbr)

hpcs.lps <- LoadH5Seurat("../data/processed/hpcs_lps_state_marked.h5Seurat")

de.state.mks <- read.csv("../outs/DE_state_markers_MAST.csv")

# C5_t2g_MF <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "MF") |>
#   dplyr::select(gs_name, entrez_gene)
# C5_t2g_BP <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP") |>
#   dplyr::select(gs_name, entrez_gene)
# C5_t2g_CC <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "CC") |>
#   dplyr::select(gs_name, entrez_gene)

# Reviewer 2
## Q2, a
# intersect(
#   intersect(
#     intersect(
#       intersect(
#         intersect(
#           de.state.mks |> filter(cell_type=="Som") |> filter(state=="Inflammation") |> pull("gene"), 
#           de.state.mks |> filter(cell_type=="Cort") |> filter(state=="Inflammation") |> pull("gene")
#         ),
#         de.state.mks |> filter(cell_type=="Lac") |> filter(state=="Inflammation") |> pull("gene")
#       ),
#       de.state.mks |> filter(cell_type=="Gonad") |> filter(state=="Inflammation") |> pull("gene")
#     ),
#     de.state.mks |> filter(cell_type=="Thyro") |> filter(state=="Inflammation") |> pull("gene")
#   ),
#   de.state.mks |> filter(cell_type=="Mel") |> filter(state=="Inflammation") |> pull("gene")
# ) -> x

# ego.bp <- groupGO(
#   gene = mapIds(
#     x = org.Mm.eg.db,
#     keys = x,
#     column = "ENTREZID",
#     keytype = "SYMBOL",
#     multiVals = "first"
#   ),
#   OrgDb = org.Mm.eg.db,
#   keyType = "ENTREZID",
#   ont = "BP",
#   level = 2,
#   readable = T
# )
# 
# ego.mf <- groupGO(
#   gene = mapIds(
#     x = org.Mm.eg.db,
#     keys = x,
#     column = "ENTREZID",
#     keytype = "SYMBOL",
#     multiVals = "first"
#   ),
#   OrgDb = org.Mm.eg.db,
#   keyType = "ENTREZID",
#   ont = "MF",
#   level = 2,
#   readable = T
# )

# Q2, b, c
cells <- LoadH5Seurat("../data/processed/cells_postprocessed.h5Seurat")
hpcs.tnf <- subset(cells, subset = treat %in% c("Saline","TNFalpha"))
hpcs.polyic <- subset(cells, subset = treat %in% c("Saline","Poly(i:c)"))

query.hpcs <- subset(cells, subset = treat %in% c("Poly(i:c)","TNFalpha")) |>
  subset(subset = cell_type_brief %in% c("Som","Lac","Cort","Gonad","Mel","Thyro")) |>
  NormalizeData() |> FindVariableFeatures() |> ScaleData() |> RunPCA()

hpcs.lps$id = "Reference"
query.hpcs$id = "Query"

som.lps <- subset(hpcs.lps, subset = cell_type_brief == "Som")
som.lps$state <- factor(som.lps$state, levels = c("Healthy","Inflammation"))
lac.lps <- subset(hpcs.lps, subset = cell_type_brief == "Lac")
lac.lps$state <- factor(lac.lps$state, levels = c("Healthy","Inflammation"))
cort.lps <- subset(hpcs.lps, subset = cell_type_brief == "Cort")
cort.lps$state <- factor(cort.lps$state, levels = c("Healthy","Inflammation"))
mel.lps <- subset(hpcs.lps, subset = cell_type_brief == "Mel")
mel.lps$state <- factor(mel.lps$state, levels = c("Healthy","Inflammation"))
gonad.lps <- subset(hpcs.lps, subset = cell_type_brief == "Gonad")
gonad.lps$state <- factor(gonad.lps$state, levels = c("Healthy","Inflammation"))
thyro.lps <- subset(hpcs.lps, subset = cell_type_brief == "Thyro")
thyro.lps$state <- factor(thyro.lps$state, levels = c("Healthy","Inflammation"))


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
query.hpcs$state <- query.hpcs$predicted.state
hpcs.merged <- merge(hpcs.lps, query.hpcs)

PADJ_CUTOFF = 0.001
FC_CUTOFF = 0.65
test_method = "MAST"

## Poly(i:c) DEGs
som.polyic <- query.hpcs |> 
  subset(subset = cell_type_brief == "Som") |>
  subset(subset = treat == "Poly(i:c)") |>
  subset(subset = predicted.state.score >= 0.9) |>
  merge(som.lps |> subset(subset = treat == "Saline"))
som.polyic$state <- factor(som.polyic$state, levels = c("Healthy","Inflammation"))
Idents(som.polyic) <- "state"
som.polyic.degs <- FindMarkers(som.polyic, ident.1 = "Inflammation", ident.2 = "Healthy", test.use = "MAST") |> 
  # dplyr::filter(p_val_adj <= PADJ_CUTOFF) %>%
  # dplyr::filter(avg_log2FC >= FC_CUTOFF | avg_log2FC <= -FC_CUTOFF) %>%
  # mutate(state = if_else(avg_log2FC > 0, "Inflammation", "Healthy")) %>%
  # relocate(state) %>%
  add_column(.before = 1, cell_type = "Som") %>%
  rownames_to_column("gene")
# entrez_ids <- mapIds(x = org.Mm.eg.db, keys = som.polyic.degs |> pull("gene"), column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
# fcs <- som.polyic.degs |> pull("avg_log2FC")
# names(fcs) = entrez_ids
# geneList = sort(fcs[!is.na(names(fcs))], decreasing = T)
# gsea_go_bp <- GSEA(geneList, TERM2GENE = C5_t2g_BP)
# gsea_go_mf <- GSEA(geneList, TERM2GENE = C5_t2g_MF)
# gsea_go_cc <- GSEA(geneList, TERM2GENE = C5_t2g_CC)

cort.polyic <- query.hpcs |> 
  subset(subset = cell_type_brief == "Cort") |>
  subset(subset = treat == "Poly(i:c)") |>
  subset(subset = predicted.state.score >= 0.9) |>
  merge(cort.lps |> subset(subset = treat == "Saline"))
cort.polyic$state <- factor(cort.polyic$state, levels = c("Healthy","Inflammation"))
Idents(cort.polyic) <- "state"
cort.polyic.degs <- FindMarkers(cort.polyic, ident.1 = "Inflammation", ident.2 = "Healthy", test.use = "MAST") |>
  dplyr::filter(p_val_adj <= PADJ_CUTOFF) %>%
  dplyr::filter(avg_log2FC >= FC_CUTOFF | avg_log2FC <= -FC_CUTOFF) %>%
  mutate(state = if_else(avg_log2FC > 0, "Inflammation", "Healthy")) %>%
  relocate(state) %>%
  add_column(.before = 1, cell_type = "Cort") %>%
  rownames_to_column("gene")

lac.polyic <- query.hpcs |> 
  subset(subset = cell_type_brief == "Lac") |>
  subset(subset = treat == "Poly(i:c)") |>
  subset(subset = predicted.state.score >= 0.9) |>
  merge(lac.lps |> subset(subset = treat == "Saline"))
lac.polyic$state <- factor(lac.polyic$state, levels = c("Healthy","Inflammation"))
Idents(lac.polyic) <- "state"
lac.polyic.degs <- FindMarkers(lac.polyic, ident.1 = "Inflammation", ident.2 = "Healthy", test.use = "MAST") |>
  dplyr::filter(p_val_adj <= PADJ_CUTOFF) %>%
  dplyr::filter(avg_log2FC >= FC_CUTOFF | avg_log2FC <= -FC_CUTOFF) %>%
  mutate(state = if_else(avg_log2FC > 0, "Inflammation", "Healthy")) %>%
  relocate(state) %>%
  add_column(.before = 1, cell_type = "Lac") %>%
  rownames_to_column("gene")

gonad.polyic <- query.hpcs |> 
  subset(subset = cell_type_brief == "Gonad") |>
  subset(subset = treat == "Poly(i:c)") |>
  subset(subset = predicted.state.score >= 0.9) |>
  merge(gonad.lps |> subset(subset = treat == "Saline"))
gonad.polyic$state <- factor(gonad.polyic$state, levels = c("Healthy","Inflammation"))
Idents(gonad.polyic) <- "state"
gonad.polyic.degs <- FindMarkers(gonad.polyic, ident.1 = "Inflammation", ident.2 = "Healthy", test.use = "MAST") |>
  dplyr::filter(p_val_adj <= PADJ_CUTOFF) %>%
  dplyr::filter(avg_log2FC >= FC_CUTOFF | avg_log2FC <= -FC_CUTOFF) %>%
  mutate(state = if_else(avg_log2FC > 0, "Inflammation", "Healthy")) %>%
  relocate(state) %>%
  add_column(.before = 1, cell_type = "Gonad") %>%
  rownames_to_column("gene")

thyro.polyic <- query.hpcs |> 
  subset(subset = cell_type_brief == "Thyro") |>
  subset(subset = treat == "Poly(i:c)") |>
  subset(subset = predicted.state.score >= 0.9) |>
  merge(thyro.lps |> subset(subset = treat == "Saline"))
thyro.polyic$state <- factor(thyro.polyic$state, levels = c("Healthy","Inflammation"))
Idents(thyro.polyic) <- "state"
thyro.polyic.degs <- FindMarkers(thyro.polyic, ident.1 = "Inflammation", ident.2 = "Healthy", test.use = "MAST") |>
  dplyr::filter(p_val_adj <= PADJ_CUTOFF) %>%
  dplyr::filter(avg_log2FC >= FC_CUTOFF | avg_log2FC <= -FC_CUTOFF) %>%
  mutate(state = if_else(avg_log2FC > 0, "Inflammation", "Healthy")) %>%
  relocate(state) %>%
  add_column(.before = 1, cell_type = "Thyro") %>%
  rownames_to_column("gene")

mel.polyic <- query.hpcs |> 
  subset(subset = cell_type_brief == "Mel") |>
  subset(subset = treat == "Poly(i:c)") |>
  subset(subset = predicted.state.score >= 0.9) |>
  merge(mel.lps |> subset(subset = treat == "Saline"))
mel.polyic$state <- factor(mel.polyic$state, levels = c("Healthy","Inflammation"))
Idents(mel.polyic) <- "state"
mel.polyic.degs <- FindMarkers(mel.polyic, ident.1 = "Inflammation", ident.2 = "Healthy", test.use = "MAST") |>
  dplyr::filter(p_val_adj <= PADJ_CUTOFF) %>%
  dplyr::filter(avg_log2FC >= FC_CUTOFF | avg_log2FC <= -FC_CUTOFF) %>%
  mutate(state = if_else(avg_log2FC > 0, "Inflammation", "Healthy")) %>%
  relocate(state) %>%
  add_column(.before = 1, cell_type = "Mel") %>%
  rownames_to_column("gene")



## TNFalpha DEGs
som.tnf <- query.hpcs |> 
  subset(subset = cell_type_brief == "Som") |>
  subset(subset = treat == "TNFalpha") |>
  subset(subset = predicted.state.score >= 0.9) |>
  merge(som.lps |> subset(subset = treat == "Saline"))
som.tnf$state <- factor(som.tnf$state, levels = c("Healthy","Inflammation"))
Idents(som.tnf) <- "state"
som.tnf.degs <- FindMarkers(som.tnf, ident.1 = "Inflammation", ident.2 = "Healthy", test.use = "MAST") |>
  # dplyr::filter(p_val_adj <= PADJ_CUTOFF) %>%
  # dplyr::filter(avg_log2FC >= FC_CUTOFF | avg_log2FC <= -FC_CUTOFF) %>%
  # mutate(state = if_else(avg_log2FC > 0, "Inflammation", "Healthy")) %>%
  # relocate(state) %>%
  add_column(.before = 1, cell_type = "Som") %>%
  rownames_to_column("gene")

cort.tnf <- query.hpcs |> 
  subset(subset = cell_type_brief == "Cort") |>
  subset(subset = treat == "TNFalpha") |>
  subset(subset = predicted.state.score >= 0.9) |>
  merge(cort.lps |> subset(subset = treat == "Saline"))
cort.tnf$state <- factor(cort.tnf$state, levels = c("Healthy","Inflammation"))
Idents(cort.tnf) <- "state"
cort.tnf.degs <- FindMarkers(cort.tnf, ident.1 = "Inflammation", ident.2 = "Healthy", test.use = "MAST")

lac.tnf <- query.hpcs |> 
  subset(subset = cell_type_brief == "Lac") |>
  subset(subset = treat == "TNFalpha") |>
  subset(subset = predicted.state.score >= 0.9) |>
  merge(lac.lps |> subset(subset = treat == "Saline"))
lac.tnf$state <- factor(lac.tnf$state, levels = c("Healthy","Inflammation"))
Idents(lac.tnf) <- "state"
lac.tnf.degs <- FindMarkers(lac.tnf, ident.1 = "Inflammation", ident.2 = "Healthy", test.use = "MAST")

gonad.tnf <- query.hpcs |> 
  subset(subset = cell_type_brief == "Gonad") |>
  subset(subset = treat == "TNFalpha") |>
  subset(subset = predicted.state.score >= 0.9) |>
  merge(gonad.lps |> subset(subset = treat == "Saline"))
gonad.tnf$state <- factor(gonad.tnf$state, levels = c("Healthy","Inflammation"))
Idents(gonad.tnf) <- "state"
gonad.tnf.degs <- FindMarkers(gonad.tnf, ident.1 = "Inflammation", ident.2 = "Healthy", test.use = "MAST")

thyro.tnf <- query.hpcs |> 
  subset(subset = cell_type_brief == "Thyro") |>
  subset(subset = treat == "TNFalpha") |>
  subset(subset = predicted.state.score >= 0.9) |>
  merge(thyro.lps |> subset(subset = treat == "Saline"))
thyro.tnf$state <- factor(thyro.tnf$state, levels = c("Healthy","Inflammation"))
Idents(thyro.tnf) <- "state"
thyro.tnf.degs <- FindMarkers(thyro.tnf, ident.1 = "Inflammation", ident.2 = "Healthy", test.use = "MAST")

mel.tnf <- query.hpcs |> 
  subset(subset = cell_type_brief == "Mel") |>
  subset(subset = treat == "TNFalpha") |>
  subset(subset = predicted.state.score >= 0.9) |>
  merge(hpcs.lps |> subset(subset = treat == "Saline"))
mel.tnf$state <- factor(mel.tnf$state, levels = c("Healthy","Inflammation"))
Idents(mel.tnf) <- "state"
mel.tnf.degs <- FindMarkers(mel.tnf, ident.1 = "Inflammation", ident.2 = "Healthy", test.use = "MAST")


library(RColorBrewer)
l.celltype.treat <- list()
l.celltype.treat[["LPS"]][["Som"]] <- som.lps; l.celltype.treat[["LPS"]][["Lac"]] <- lac.lps; l.celltype.treat[["LPS"]][["Cort"]] <- cort.lps
l.celltype.treat[["LPS"]][["Mel"]] <- mel.lps; l.celltype.treat[["LPS"]][["Gonad"]] <- gonad.lps; l.celltype.treat[["LPS"]][["Thyro"]] <- thyro.lps
l.celltype.treat[["Poly(i:c)"]][["Som"]] <- som.polyic; l.celltype.treat[["Poly(i:c)"]][["Lac"]] <- lac.polyic; l.celltype.treat[["Poly(i:c)"]][["Cort"]] <- cort.polyic
l.celltype.treat[["Poly(i:c)"]][["Mel"]] <- mel.polyic; l.celltype.treat[["Poly(i:c)"]][["Gonad"]] <- gonad.polyic; l.celltype.treat[["Poly(i:c)"]][["Thyro"]] <- thyro.polyic
l.celltype.treat[["TNF"]][["Som"]] <- som.tnf; l.celltype.treat[["TNF"]][["Lac"]] <- lac.tnf; l.celltype.treat[["TNF"]][["Cort"]] <- cort.tnf
l.celltype.treat[["TNF"]][["Mel"]] <- mel.tnf; l.celltype.treat[["TNF"]][["Gonad"]] <- gonad.tnf; l.celltype.treat[["TNF"]][["Thyro"]] <- thyro.tnf


for(){
  
}

VlnPlot(som.lps, features = "Serpina3n", group.by = "state", cols = brewer.pal(n = 12, name = "Paired")[c(2,8)]) +
  theme(axis.title = element_blank(), legend.position = "none", plot.title = element_blank(), axis.text.x = element_blank())
ggsave("Serpina3n_som_lps.eps", device = "eps", path = "../figures/FigS2/violinplots/", width = 4, height = 3, dpi = 300)
VlnPlot(lac.lps, features = "Serpina3n", group.by = "state", cols = brewer.pal(n = 12, name = "Paired")[c(2,8)]) +
  theme(axis.title = element_blank(), legend.position = "none", plot.title = element_blank(), axis.text.x = element_blank())
ggsave("Serpina3n_lac_lps.eps", device = "eps", path = "../figures/FigS2/violinplots/", width = 4, height = 3, dpi = 300)
VlnPlot(cort.lps, features = "Serpina3n", group.by = "state", cols = brewer.pal(n = 12, name = "Paired")[c(2,8)]) +
  theme(axis.title = element_blank(), legend.position = "none", plot.title = element_blank(), axis.text.x = element_blank())
ggsave("Serpina3n_cort_lps.eps", device = "eps", path = "../figures/FigS2/violinplots/", width = 4, height = 3, dpi = 300)
VlnPlot(gonad.lps, features = "Serpina3n", group.by = "state", cols = brewer.pal(n = 12, name = "Paired")[c(2,8)]) +
  theme(axis.title = element_blank(), legend.position = "none", plot.title = element_blank(), axis.text.x = element_blank())
ggsave("Serpina3n_gonad_lps.eps", device = "eps", path = "../figures/FigS2/violinplots/", width = 4, height = 3, dpi = 300)
VlnPlot(mel.lps, features = "Serpina3n", group.by = "state", cols = brewer.pal(n = 12, name = "Paired")[c(2,8)]) +
  theme(axis.title = element_blank(), legend.position = "none", plot.title = element_blank(), axis.text.x = element_blank())
ggsave("Serpina3n_mel_lps.eps", device = "eps", path = "../figures/FigS2/violinplots/", width = 4, height = 3, dpi = 300)
VlnPlot(thyro.lps, features = "Serpina3n", group.by = "state", cols = brewer.pal(n = 12, name = "Paired")[c(2,8)]) +
  theme(axis.title = element_blank(), legend.position = "none", plot.title = element_blank(), axis.text.x = element_blank())
ggsave("Serpina3n_thyro_lps.eps", device = "eps", path = "../figures/FigS2/violinplots/", width = 4, height = 3, dpi = 300)

VlnPlot(som.polyic, features = "Serpina3n", group.by = "state", cols = brewer.pal(n = 12, name = "Paired")[c(2,8)]) +
  theme(axis.title = element_blank(), legend.position = "none", plot.title = element_blank(), axis.text.x = element_blank())
ggsave("Serpina3n_som_lps.eps", device = "eps", path = "../figures/FigS2/violinplots/", width = 4, height = 3, dpi = 300)
VlnPlot(lac.lps, features = "Serpina3n", group.by = "state", cols = brewer.pal(n = 12, name = "Paired")[c(2,8)]) +
  theme(axis.title = element_blank(), legend.position = "none", plot.title = element_blank(), axis.text.x = element_blank())
ggsave("Serpina3n_lac_lps.eps", device = "eps", path = "../figures/FigS2/violinplots/", width = 4, height = 3, dpi = 300)
VlnPlot(cort.lps, features = "Serpina3n", group.by = "state", cols = brewer.pal(n = 12, name = "Paired")[c(2,8)]) +
  theme(axis.title = element_blank(), legend.position = "none", plot.title = element_blank(), axis.text.x = element_blank())
ggsave("Serpina3n_cort_lps.eps", device = "eps", path = "../figures/FigS2/violinplots/", width = 4, height = 3, dpi = 300)
VlnPlot(gonad.lps, features = "Serpina3n", group.by = "state", cols = brewer.pal(n = 12, name = "Paired")[c(2,8)]) +
  theme(axis.title = element_blank(), legend.position = "none", plot.title = element_blank(), axis.text.x = element_blank())
ggsave("Serpina3n_gonad_lps.eps", device = "eps", path = "../figures/FigS2/violinplots/", width = 4, height = 3, dpi = 300)
VlnPlot(mel.lps, features = "Serpina3n", group.by = "state", cols = brewer.pal(n = 12, name = "Paired")[c(2,8)]) +
  theme(axis.title = element_blank(), legend.position = "none", plot.title = element_blank(), axis.text.x = element_blank())
ggsave("Serpina3n_mel_lps.eps", device = "eps", path = "../figures/FigS2/violinplots/", width = 4, height = 3, dpi = 300)
VlnPlot(thyro.lps, features = "Serpina3n", group.by = "state", cols = brewer.pal(n = 12, name = "Paired")[c(2,8)]) +
  theme(axis.title = element_blank(), legend.position = "none", plot.title = element_blank(), axis.text.x = element_blank())
ggsave("Serpina3n_thyro_lps.eps", device = "eps", path = "../figures/FigS2/violinplots/", width = 4, height = 3, dpi = 300)



VlnPlot(hpcs.polyic, features = "Serpina3n", group.by = "state", cols = brewer.pal(n = 12, name = "Paired")[c(2,8)]) +
  theme(axis.title = element_blank(), legend.position = "none", plot.title = element_blank())
VlnPlot(hpcs.tnf, features = "Serpina3n", group.by = "state", cols = brewer.pal(n = 12, name = "Paired")[c(2,8)]) +
  theme(axis.title = element_blank(), legend.position = "none", plot.title = element_blank())






# library(plot1cell)
# hpcs.merged$cell_type_brief <- factor(hpcs.merged$cell_type_brief, levels = c("Som", "Lac", "Cort", "Mel", "Gonad", "Thyro"))
# Idents(hpcs.merged) <- "cell_type_brief"

# complex_dotplot_multiple(
#   hpcs.lps, features = c("Cxcl1","Cxcl10","Serpina3n","Irf7"), groups = "state", celltypes = c("Som", "Lac", "Cort", "Mel", "Gonad", "Thyro")
# )
# 
# complex_dotplot_multiple(
#   hpcs.polyic, features = c("Cxcl1","Cxcl10","Serpina3n","Irf7"), groups = "state", celltypes = c("Som", "Lac", "Cort", "Mel", "Gonad", "Thyro")
# )
# 
# complex_dotplot_multiple(
#   hpcs.tnf, features = c("Cxcl1","Cxcl10","Serpina3n","Irf7"), groups = "state", celltypes = c("Som", "Lac", "Cort", "Mel", "Gonad", "Thyro")
# )

## Q3, d
VlnPlot(cort.polyic, features = "Nptx2", group.by = "state", cols = brewer.pal(n = 12, name = "Paired")[c(2,8)]) +
  theme(axis.title = element_blank(), legend.position = "none", plot.title = element_blank(), axis.text.x = element_blank())
ggsave("Nptx2_cort_polyic.eps", device = "eps", path = "../figures/FigS2/violinplots/", width = 4, height = 3, dpi = 300)
VlnPlot(cort.tnf, features = "Nptx2", group.by = "state", cols = brewer.pal(n = 12, name = "Paired")[c(2,8)]) +
  theme(axis.title = element_blank(), legend.position = "none", plot.title = element_blank(), axis.text.x = element_blank())
ggsave("Nptx2_cort_tnf.eps", device = "eps", path = "../figures/FigS2/violinplots/", width = 4, height = 3, dpi = 300)


