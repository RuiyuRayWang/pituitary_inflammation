library(Seurat)
library(SeuratDisk)
library(ggpubr)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

suppressMessages(
  extrafont::loadfonts(device="postscript")
)

hpcs.lps <- LoadH5Seurat("../data/processed/hpcs_lps_state_marked.h5Seurat")

# GO analysis
de.state.mks <- read.csv(file = '../outs/de_state_markers.csv')
de.state.mks.uniq <- de.state.mks[!duplicated(de.state.mks$gene),]

de.state.mks.uniq$entrez <- mapIds(
  x = org.Mm.eg.db, 
  keys = de.state.mks.uniq$gene, 
  column = "ENTREZID", 
  keytype = "SYMBOL", 
  multiVals = "first"
)

input_genes <- de.state.mks.uniq %>% dplyr::filter(state == "Inflammation") %>% pull(entrez) %>% na.omit()
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
    legend.text = element_text(size = 12),
    axis.text.y = element_text(size = 15)
  )

GO_CC_terms_use <- c("GO:0033646","GO:0043657","GO:0005615","GO:0043230","GO:0042824","GO:0030670","GO:0031410","GO:0042612","GO:0009897","GO:0030139","GO:0005789","GO:0048471")
de_ego.CC.filtered <- de_ego.CC %>% dplyr::filter(p.adjust < 0.05) %>% dplyr::filter(ID %in% GO_CC_terms_use)
df_tmp <- de_ego.CC.filtered@result

df_tmp <- dplyr::arrange(df_tmp, DOSE::parse_ratio(GeneRatio))
y_col <- ifelse(df_tmp$ID %in% c("GO:0005615","GO:0005576","GO:0099503"), "red", "black")

p2 <- clusterProfiler::dotplot(de_ego.CC.filtered, showCategory = nrow(de_ego.CC.filtered)) +
  theme(
    legend.text = element_text(size = 12),
    axis.text.y=element_text(color=y_col, size = 15)
  )

p <- ggarrange(p1, p2, ncol = 1, nrow = 2, common.legend = TRUE, legend = "right")

ggsave(
  filename = "GO.dotplot.eps", 
  plot = p,
  device = "eps", 
  path = '../figures/Fig2/', 
  width = 180, height = 250, 
  dpi = 300, 
  units = "mm",
  family = "Arial"
)