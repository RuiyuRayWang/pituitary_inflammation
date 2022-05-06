library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

suppressMessages(
  extrafont::loadfonts(device="postscript")
)

hpcs.lps <- LoadH5Seurat("data/hpcs_lps_state_marked.h5Seurat")

## Conserved Markers
MIN_P_CUTOFF = 0.01
FC_CUTOFF = 0.8

Idents(hpcs.lps) <- "cell_type_brief"
som.csvd.mks <- FindConservedMarkers(object = hpcs.lps, ident.1 = "Som", grouping.var = "state", only.pos = TRUE) %>%
  dplyr::filter(minimump_p_val <= MIN_P_CUTOFF) %>%
  dplyr::filter(Healthy_avg_log2FC >= FC_CUTOFF) %>% 
  dplyr::filter(Inflammation_avg_log2FC >= FC_CUTOFF) %>%
  add_column(.before = 1, cell_type = "Som") %>% 
  rownames_to_column("gene")
lac.csvd.mks <- FindConservedMarkers(object = hpcs.lps, ident.1 = "Lac", grouping.var = "state", only.pos = TRUE) %>%
  dplyr::filter(minimump_p_val <= MIN_P_CUTOFF) %>%
  dplyr::filter(Healthy_avg_log2FC >= FC_CUTOFF) %>% 
  dplyr::filter(Inflammation_avg_log2FC >= FC_CUTOFF) %>%
  add_column(.before = 1, cell_type = "Lac") %>% 
  rownames_to_column("gene")
cort.csvd.mks <- FindConservedMarkers(object = hpcs.lps, ident.1 = "Cort", grouping.var = "state", only.pos = TRUE) %>%
  dplyr::filter(minimump_p_val <= MIN_P_CUTOFF) %>%
  dplyr::filter(Healthy_avg_log2FC >= FC_CUTOFF) %>% 
  dplyr::filter(Inflammation_avg_log2FC >= FC_CUTOFF) %>%
  add_column(.before = 1, cell_type = "Cort") %>% 
  rownames_to_column("gene")
mel.csvd.mks <- FindConservedMarkers(object = hpcs.lps, ident.1 = "Mel", grouping.var = "state", only.pos = TRUE) %>%
  dplyr::filter(minimump_p_val <= MIN_P_CUTOFF) %>%
  dplyr::filter(Healthy_avg_log2FC >= FC_CUTOFF) %>% 
  dplyr::filter(Inflammation_avg_log2FC >= FC_CUTOFF) %>%
  add_column(.before = 1, cell_type = "Mel") %>% 
  rownames_to_column("gene")
gonad.csvd.mks <- FindConservedMarkers(object = hpcs.lps, ident.1 = "Gonad", grouping.var = "state", only.pos = TRUE) %>%
  dplyr::filter(minimump_p_val <= MIN_P_CUTOFF) %>%
  dplyr::filter(Healthy_avg_log2FC >= FC_CUTOFF) %>% 
  dplyr::filter(Inflammation_avg_log2FC >= FC_CUTOFF) %>%
  add_column(.before = 1, cell_type = "Gonad") %>% 
  rownames_to_column("gene")
thyro.csvd.mks <- FindConservedMarkers(object = hpcs.lps, ident.1 = "Thyro", grouping.var = "state", only.pos = TRUE) %>%
  dplyr::filter(minimump_p_val <= MIN_P_CUTOFF) %>%
  dplyr::filter(Healthy_avg_log2FC >= FC_CUTOFF) %>% 
  dplyr::filter(Inflammation_avg_log2FC >= FC_CUTOFF) %>%
  add_column(.before = 1, cell_type = "Thyro") %>% 
  rownames_to_column("gene")

csvd.cell.mks <- bind_rows(som.csvd.mks, lac.csvd.mks) %>% 
  bind_rows(cort.csvd.mks) %>% 
  bind_rows(mel.csvd.mks) %>% 
  bind_rows(gonad.csvd.mks) %>%
  bind_rows(thyro.csvd.mks)

### Note there are duplicated values
csvd.cell.mks %>% write.csv(file = "../outs/conserved_cell_markers.csv", row.names = FALSE)

## Conserved inflammatory hallmark genes
MIN_P_CUTOFF = 0.01
PADJ_CUTOFF = 0.1

Idents(hpcs.lps) <- "state"
csvd.state.mks <- FindConservedMarkers(
  object = hpcs.lps, 
  ident.1 = "Inflammation", 
  grouping.var = "cell_type_brief"
)

csvd.state.mks <- csvd.state.mks %>% 
  dplyr::filter(minimump_p_val <= MIN_P_CUTOFF) %>%
  dplyr::filter(Som_p_val_adj <= PADJ_CUTOFF) %>%
  dplyr::filter(Lac_p_val_adj <= PADJ_CUTOFF) %>%
  dplyr::filter(Cort_p_val_adj <= PADJ_CUTOFF) %>%
  dplyr::filter(Mel_p_val_adj <= PADJ_CUTOFF) %>%
  dplyr::filter(Gonad_p_val_adj <= PADJ_CUTOFF) %>%
  dplyr::filter(Thyro_p_val_adj <= PADJ_CUTOFF) %>%
  rownames_to_column("gene")

csvd.state.mks <- csvd.state.mks %>% 
  mutate(state = 
           if_else(
             Som_avg_log2FC>0 & Lac_avg_log2FC>0 & Cort_avg_log2FC>0 & Mel_avg_log2FC>0 & Gonad_avg_log2FC>0 & Thyro_avg_log2FC>0,
             "Inflammation", 
             if_else(Som_avg_log2FC<0 & Lac_avg_log2FC<0 & Cort_avg_log2FC<0 & Mel_avg_log2FC<0 & Gonad_avg_log2FC<0 & Thyro_avg_log2FC<0,
                     "Healthy",
                     "Ambiguous")
             )
         ) %>%
  relocate(state)
csvd.state.mks %>% write.csv(file = "../outs/conserved_state_markers.csv", row.names = FALSE)

## Differentially Expressed
PADJ_CUTOFF = 0.01
FC_CUTOFF = 0.8

Idents(hpcs.lps) <- "state"
som.state.mks <- subset(hpcs.lps, subset = cell_type_brief == "Som") %>% FindMarkers(ident.1 = "Inflammation") %>%
  dplyr::filter(p_val_adj <= PADJ_CUTOFF) %>%
  dplyr::filter(avg_log2FC >= FC_CUTOFF | avg_log2FC <= -FC_CUTOFF) %>%
  mutate(state = if_else(avg_log2FC > 0, "Inflammation", "Healthy")) %>%
  relocate(state) %>%
  add_column(.before = 1, cell_type = "Som") %>%
  rownames_to_column("gene")
lac.state.mks <- subset(hpcs.lps, subset = cell_type_brief == "Lac") %>% FindMarkers(ident.1 = "Inflammation") %>%
  dplyr::filter(p_val_adj <= PADJ_CUTOFF) %>%
  dplyr::filter(avg_log2FC >= FC_CUTOFF | avg_log2FC <= -FC_CUTOFF) %>%
  mutate(state = if_else(avg_log2FC > 0, "Inflammation", "Healthy")) %>%
  relocate(state) %>%
  add_column(.before = 1, cell_type = "Lac") %>%
  rownames_to_column("gene")
cort.state.mks <- subset(hpcs.lps, subset = cell_type_brief == "Cort") %>% FindMarkers(ident.1 = "Inflammation") %>%
  dplyr::filter(p_val_adj <= PADJ_CUTOFF) %>%
  dplyr::filter(avg_log2FC >= FC_CUTOFF | avg_log2FC <= -FC_CUTOFF) %>%
  mutate(state = if_else(avg_log2FC > 0, "Inflammation", "Healthy")) %>%
  relocate(state) %>%
  add_column(.before = 1, cell_type = "Cort") %>%
  rownames_to_column("gene")
mel.state.mks <- subset(hpcs.lps, subset = cell_type_brief == "Mel") %>% FindMarkers(ident.1 = "Inflammation") %>%
  dplyr::filter(p_val_adj <= PADJ_CUTOFF) %>%
  dplyr::filter(avg_log2FC >= FC_CUTOFF | avg_log2FC <= -FC_CUTOFF) %>%
  mutate(state = if_else(avg_log2FC > 0, "Inflammation", "Healthy")) %>%
  relocate(state) %>%
  add_column(.before = 1, cell_type = "Mel") %>%
  rownames_to_column("gene")
gonad.state.mks <- subset(hpcs.lps, subset = cell_type_brief == "Gonad") %>% FindMarkers(ident.1 = "Inflammation") %>%
  dplyr::filter(p_val_adj <= PADJ_CUTOFF) %>%
  dplyr::filter(avg_log2FC >= FC_CUTOFF | avg_log2FC <= -FC_CUTOFF) %>%
  mutate(state = if_else(avg_log2FC > 0, "Inflammation", "Healthy")) %>%
  relocate(state) %>%
  add_column(.before = 1, cell_type = "Gonad") %>%
  rownames_to_column("gene")
thyro.state.mks <- subset(hpcs.lps, subset = cell_type_brief == "Thyro") %>% FindMarkers(ident.1 = "Inflammation") %>%
  dplyr::filter(p_val_adj <= PADJ_CUTOFF) %>%
  dplyr::filter(avg_log2FC >= FC_CUTOFF | avg_log2FC <= -FC_CUTOFF) %>%
  mutate(state = if_else(avg_log2FC > 0, "Inflammation", "Healthy")) %>%
  relocate(state) %>%
  add_column(.before = 1, cell_type = "Thyro") %>%
  rownames_to_column("gene")

de.state.mks <- bind_rows(som.state.mks, lac.state.mks) %>% 
  bind_rows(cort.state.mks) %>% 
  bind_rows(mel.state.mks) %>% 
  bind_rows(gonad.state.mks) %>%
  bind_rows(thyro.state.mks)

de.state.mks %>% write.csv(file = "../outs/de_state_markers.csv", row.names = FALSE)


# Plot Heatmap
source("DoMultiBarHeatmap.R")

## Conserved Markers
### Plot order
csvd.cell.mks$cell_type <- factor(csvd.cell.mks$cell_type, levels = c("Som","Lac","Cort","Mel","Gonad","Thyro"))

### Purge duplilcates
csvd.cell.mks <- csvd.cell.mks[!duplicated(csvd.cell.mks$gene),]

### Gene scaling
hpcs.lps <- ScaleData(hpcs.lps, features = union(VariableFeatures(hpcs.lps), csvd.cell.mks$gene))

### Make Plot
heatmap.csvd.cells <- DoMultiBarHeatmap(
  hpcs.lps,
  features = csvd.cell.mks$gene,
  label = FALSE,
  group.bar = FALSE,
  group.by = 'cell_type_brief',
  additional.group.by = 'state',
  additional.group.sort.by = 'state'
) + 
  scale_fill_distiller(palette = 'RdYlBu') + 
  NoAxes()
ggsave(
  filename = "heatmap_csvd_cells_hpcslps.eps",
  plot = heatmap.csvd.cells, 
  device = "eps", 
  path = "../figures/Fig1/", 
  dpi = 300
)
heatmap.csvd.cells.clean <- DoMultiBarHeatmap(
  hpcs.lps,
  features = csvd.cell.mks$gene,
  label = FALSE,
  group.bar = FALSE,
  group.by = 'cell_type_brief',
  additional.group.by = 'state',
  additional.group.sort.by = 'state'
  ) + 
  scale_fill_distiller(palette = 'RdYlBu') +
  NoLegend() + 
  NoAxes()  # Remove gene names
ggsave(
  filename = "heatmap_csvd_cells_hpcslps_clean.eps",
  plot = heatmap.csvd.cells.clean, 
  device = "eps", 
  path = "../figures/Fig1/", 
  dpi = 300
)

# ## Conserved DE hallmark genes between states
# ### Gene scaling
# hpcs.lps <- ScaleData(hpcs.lps, features = union(VariableFeatures(hpcs.lps), csvd.state.mks$gene))
# 
# ### Make Plots
# heatmap.csvd.state <- DoMultiBarHeatmap(
#   hpcs.lps,
#   features = csvd.state.mks$gene,
#   label = FALSE,
#   group.bar = FALSE,
#   group.by = 'cell_type_brief',
#   additional.group.by = 'state',
#   additional.group.sort.by = 'state'
# ) + 
#   scale_fill_distiller(palette = 'PRGn') + 
#   NoAxes()
# ggsave(
#   filename = "heatmap_csvd_state_hpcslps.eps",
#   plot = heatmap.csvd.state, 
#   device = "eps", 
#   path = "../figures/Fig1/", 
#   dpi = 300
# )
# heatmap.csvd.state.clean <- DoMultiBarHeatmap(
#   hpcs.lps,
#   features = csvd.state.mks$gene,
#   label = FALSE,
#   group.bar = FALSE,
#   group.by = 'cell_type_brief',
#   additional.group.by = 'state',
#   additional.group.sort.by = 'state'
# ) + 
#   scale_fill_distiller(palette = 'PRGn') +
#   NoLegend() + 
#   NoAxes()  # Remove gene names
# ggsave(
#   filename = "heatmap_csvd_state_hpcslps_clean.eps",
#   plot = heatmap.csvd.state.clean, 
#   device = "eps", 
#   path = "../figures/Fig1/", 
#   dpi = 300
# )

## Differentially Expressed Markers
de.state.mks$cell_type <- factor(de.state.mks$cell_type, levels = c("Som","Lac","Cort","Mel","Gonad","Thyro"))
de.state.mks$state <- factor(de.state.mks$state, levels = c("Healthy","Inflammation"))

### Purge duplilcates
state.mks.uniq <- de.state.mks[!duplicated(de.state.mks$gene),]

### Gene scaling
hpcs.lps <- ScaleData(hpcs.lps, features = union(VariableFeatures(hpcs.lps), state.mks.uniq$gene))

### Make Plot
heatmap.de <- DoMultiBarHeatmap(
  hpcs.lps,
  features = state.mks.uniq %>% arrange(cell_type, state) %>% pull(gene),
  label = FALSE,
  group.bar = FALSE,
  group.by = 'cell_type_brief',
  additional.group.by = 'state',
  additional.group.sort.by = 'state') + 
  scale_fill_distiller(palette = 'PiYG') +
  NoAxes()
ggsave(
  filename = "heatmap_state_hpcslps.eps",
  plot = heatmap.de, 
  device = "eps", 
  path = "../figures/Fig1/", 
  dpi = 300
)
heatmap.de.clean <- DoMultiBarHeatmap(
  hpcs.lps,
  features = state.mks.uniq %>% arrange(cell_type, state) %>% pull(gene),
  label = FALSE,
  group.bar = FALSE,
  group.by = 'cell_type_brief',
  additional.group.by = 'state',
  additional.group.sort.by = 'state') + 
  scale_fill_distiller(palette = 'PiYG') +
  NoLegend() + NoAxes()  # Clean plot
ggsave(
  filename = "heatmap_state_hpcslps_clean.eps",
  plot = heatmap.de.clean, 
  device = "eps", 
  path = "../figures/Fig1/", 
  dpi = 300
)

# GO analysis
de.state.mks <- read.csv(file = '../outs/de_state_markers.csv')
de.state.mks.uniq <- de.state.mks[!duplicated(de.state.mks$gene),]

state.mks.uniq$entrez <- mapIds(
  x = org.Mm.eg.db, 
  keys = state.mks.uniq$gene, 
  column = "ENTREZID", 
  keytype = "SYMBOL", 
  multiVals = "first"
)

input_genes <- state.mks.uniq %>% dplyr::filter(state == "Inflammation") %>% pull(entrez) %>% na.omit()
background <- mapIds(
  x = org.Mm.eg.db, 
  keys = rownames(hpcs.lps), 
  column = "ENTREZID", 
  keytype = "SYMBOL", 
  multiVals = "first"
) %>% na.omit()

## Cellular Component (CC)
de_ego.CC <- enrichGO(
  gene     = input_genes,
  universe = background,
  OrgDb    = org.Mm.eg.db,
  ont      = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  minGSSize = 10,
  maxGSSize = 2000,  # To get [GO:0005615] extracellular space and [GO:0005576] extracellular region, tune this value
  readable = TRUE
  )
## Molecular Function (MF)
de_ego.MF <- enrichGO(
  gene     = input_genes,
  universe = background,
  OrgDb    = org.Mm.eg.db,
  ont      = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)
## Biological Process (BP)
de_ego.BP <- enrichGO(
  gene     = input_genes,
  universe = background,
  OrgDb    = org.Mm.eg.db,
  ont      = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

## Make Plots
GO_BP_terms_use <- c("GO:0035456","GO:0034341","GO:0019221","GO:0031349","GO:0002237","GO:0071222","GO:0045088","GO:0050727")
de_ego.BP.filtered <- de_ego.BP %>% dplyr::filter(p.adjust < 0.05) %>% dplyr::filter(ID %in% GO_BP_terms_use)

p1 <- clusterProfiler::dotplot(de_ego.BP.filtered, showCategory = nrow(de_ego.BP.filtered)) +
  theme(
    axis.text.y=element_text(size = 11)
  )

GO_CC_terms_use <- c("GO:0033646","GO:0043657","GO:0005615","GO:0043230","GO:0042824","GO:0030670","GO:0031410","GO:0042612","GO:0009897","GO:0030139","GO:0005789","GO:0048471")
de_ego.CC.filtered <- de_ego.CC %>% dplyr::filter(p.adjust < 0.05) %>% dplyr::filter(ID %in% GO_CC_terms_use)
df_tmp <- de_ego.CC.filtered@result

df_tmp <- dplyr::arrange(df_tmp, DOSE::parse_ratio(GeneRatio))
y_col <- ifelse(df_tmp$ID %in% c("GO:0005615","GO:0005576","GO:0099503"), "red", "black")

p2 <- clusterProfiler::dotplot(de_ego.CC.filtered, showCategory = nrow(de_ego.CC.filtered)) +
  theme(
    axis.text.y=element_text(color=y_col, size = 11)
  )

p <- ggarrange(p1, p2, common.legend = TRUE, legend = "right")

ggsave(
  filename = "GO.dotplot.eps", 
  plot = p,
  device = "eps", 
  path = '../figures/Fig1/', 
  width = 240, height = 90, 
  dpi = 300, 
  units = "mm", 
  family = "Arial"
  )
