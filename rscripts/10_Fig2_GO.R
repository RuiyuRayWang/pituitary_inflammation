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

test_method = "MAST"

# GO analysis
de.state.mks <- read.csv(file = paste0("../outs/de_state_markers_",test_method,".csv"))
de.state.mks.uniq <- de.state.mks[!duplicated(de.state.mks$gene),]
# de.state.mks$initial <- substr(de.state.mks$cell_type, start=1, stop =1)
# de.state.mks <- de.state.mks %>% dplyr::group_by(gene) %>% dplyr::mutate(venn_relation = paste0(initial, collapse = ":")) %>% dplyr::select(-c(initial)) %>% dplyr::distinct()

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
de_ego.CC@result %>% write.csv("../outs/GO_CC.csv", row.names = FALSE, quote = TRUE)
saveRDS(de_ego.CC, "../data/processed/GO_CC.rds")

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
de_ego.MF@result %>% write.csv("../outs/GO_MF.csv", row.names = FALSE, quote = TRUE)
saveRDS(de_ego.MF, "../data/processed/GO_MF.rds")

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
de_ego.BP@result %>% write.csv("../outs/GO_BP.csv", row.names = FALSE, quote = TRUE)
saveRDS(de_ego.BP, "../data/processed/GO_BP.rds")

## Make Plots
de_ego.BP <- readRDS("../data/processed/GO_BP.rds")
de_ego.CC <- readRDS("../data/processed/GO_CC.rds")
de_ego.MF <- readRDS("../data/processed/GO_MF.rds")
### BP
GO_BP_terms_use <- c("GO:0035456","GO:0034341","GO:0019221","GO:0031349","GO:0002237","GO:0071222","GO:0045088","GO:0050727")
de_ego.BP.filtered <- de_ego.BP %>% dplyr::filter(p.adjust < 0.05) %>% dplyr::filter(ID %in% GO_BP_terms_use)

p1 <- clusterProfiler::dotplot(de_ego.BP.filtered, showCategory = nrow(de_ego.BP.filtered)) +
  scale_x_continuous(
    breaks = c(0.06,0.07,0.08,0.09)
      ) +
  scale_y_discrete(
    labels=function(x) str_wrap(x, width=20)
    ) +
  scale_size(breaks = c(25,35)) +
  scale_color_continuous(
    low = "red",
    high = "blue",
    breaks = c(2e-9, 1e-10),
    guide = guide_colorbar(reverse=TRUE, nbin = 500)
  ) +
  guides(
    color = guide_colorbar(order = 1, reverse = TRUE),
    size = guide_legend(order = 2)
  ) +
  theme(
    legend.text = element_text(size = 12, hjust = .1),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 15), 
    axis.text.y = element_text(size = 16, lineheight = .8),
    legend.position = c(.1,-.18),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    plot.margin = unit(c(.2,.2,2,.5),"cm")
  )

### CC
GO_CC_terms_use <- c("GO:0033646","GO:0043657","GO:0005615","GO:0043230","GO:0042824","GO:0030670","GO:0031410","GO:0042612","GO:0009897")  ## "GO:0030139","GO:0005789","GO:0048471"
de_ego.CC.filtered <- de_ego.CC %>% dplyr::filter(p.adjust < 0.05) %>% dplyr::filter(ID %in% GO_CC_terms_use)
# df_tmp <- de_ego.CC.filtered@result
# 
# df_tmp <- dplyr::arrange(df_tmp, DOSE::parse_ratio(GeneRatio))
# y_col <- ifelse(df_tmp$ID %in% c("GO:0005615","GO:0005576","GO:0099503"), "red", "black")

p2 <- clusterProfiler::dotplot(de_ego.CC.filtered, showCategory = nrow(de_ego.CC.filtered)) +
  scale_y_discrete(
    labels=function(x) str_wrap(x, width=20)
    ) +
  scale_size(breaks = c(20,60)) +
  scale_color_continuous(
    low = "red",
    high = "blue",
    breaks = c(0.002, 0.0002),
    guide = guide_colorbar(reverse=TRUE, nbin = 500)
  ) +
  guides(
    color = guide_colorbar(order = 1, reverse = TRUE),
    size = guide_legend(order = 2)
  ) +
  theme(
    legend.text = element_text(size = 12, hjust = .1),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 15), 
    axis.text.y = element_text(size = 16, lineheight = .8),
    legend.position = c(.1,-.18),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    plot.margin = unit(c(.2,.2,2,.5),"cm")
  )

ggsave(
  filename = "GO_BP_dotplot.eps",
  plot = p1,
  device = "eps",
  path = '../figures/Fig2/',
  width = 120, height = 150,
  dpi = 300,
  units = "mm",
  family = "Arial"
)

ggsave(
  filename = "GO_CC_dotplot.eps",
  plot = p2,
  device = "eps",
  path = '../figures/Fig2/',
  width = 120, height = 150,
  dpi = 300,
  units = "mm",
  family = "Arial"
)

# DEPRECATED
# p <- ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")  ## BUG: inconsistent p.adjust color for the shared legend
# p <- p1 | p2

# ggsave(
#   filename = "GO.dotplot.eps", 
#   plot = p,
#   device = "eps", 
#   path = '../figures/Fig2/', 
#   width = 280, height = 150, 
#   dpi = 300, 
#   units = "mm",
#   family = "Arial"
# )