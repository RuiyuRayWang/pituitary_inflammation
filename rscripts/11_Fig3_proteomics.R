library(clusterProfiler)
library(org.Mm.eg.db)
library(EnhancedVolcano)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

suppressMessages(
  extrafont::loadfonts(device="postscript")
)

ms_dat <- read.csv(file = "../data/proteomics/ms_dat.csv")

# Customized coloring scheme: 
#   4.9 Over-ride colouring scheme with custom key-value pairs
# <https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html#over-ride-colouring-scheme-with-custom-key-value-pairs>

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
  ms_dat$FC < -1 & ms_dat$Pvalue < 10^-1.3, 'royalblue',
  ifelse(ms_dat$FC > 1 & ms_dat$Pvalue < 10^-1.3, 'orange',
         'grey'))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'orange'] <- 'up'
names(keyvals)[keyvals == 'grey'] <- 'no'
names(keyvals)[keyvals == 'royalblue'] <- 'down'

FC_cutoff = 1
p_cutoff = 10^-1.3
EnhancedVolcano(ms_dat,
                lab = ms_dat$gene,
                selectLab = c("Hmox1","Vegfa","Eif2ak1","Eif2ak3","Dusp1","Trib3","Ddit4","Ddit3","Ppp1r15a"),
                x = "FC",
                y = "Pvalue",
                xlim = c(-8,9),
                ylim = c(0,7),
                labSize = 4,
                pointSize = 2.5,
                pCutoff = 10^-1.3,
                FCcutoff = 1,
                colCustom = keyvals,
                colAlpha = 1,
                xlab = bquote(~Log[2]~ 'fold change: LPS/Saline'),
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                drawConnectors = TRUE,
                widthConnectors = 0.6,
                arrowhead = FALSE,
                colConnectors = "black",
                legendPosition = "right"
) +
  theme(panel.border = element_rect(fill = NA))

ggsave(
  filename = "volcano_plot_proteomics.eps", 
  device = "eps", 
  path = "../figures/Fig3/", 
  width = 160, height = 120, 
  dpi = 300, 
  units = "mm", 
  family = "Arial"
  )

ms_ego.CC <- enrichGO(gene     = mapIds(org.Mm.eg.db,
                                        keys = ms_dat[ms_dat$sig=="up","ensembl_id"],
                                        column = "ENTREZID",
                                        keytype = "ENSEMBL",
                                        multiVals = "first"),
                      universe = as.character(mapIds(org.Mm.eg.db,
                                                     keys = ms_dat$ensembl_id,
                                                     column = "ENTREZID",
                                                     keytype = "ENSEMBL",
                                                     multiVals = "first")),
                      OrgDb    = org.Mm.eg.db,
                      ont      = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.4,
                      qvalueCutoff = 0.4,
                      readable = TRUE
)


y_col <- rev(c("red","red","black","black","black"))

clusterProfiler::dotplot(ms_ego.CC, showCategory = 5, color = "pvalue") + 
  scale_color_gradient(low = "#CCCC00", high = "#FF6600", trans = "reverse") +
  theme(
    axis.text.y = element_text(color=y_col),
    axis.text.x = element_text(size = 10)
    )

ggsave(
  filename = "GO_CC_dotplot_proteomics.eps", 
  device = "eps", 
  path = "../figures/Fig3/", 
  width = 120, height = 90, 
  dpi = 300, 
  units = "mm", 
  family = "Arial"
  )

ms.GO_CC.df <- ms_ego.CC@result %>% slice_head(n = 5)
# barplot(ms_ego.CC, showCategory = 5, color = "pvalue", order = T) + 
#   scale_fill_gradient(low = "#CCCC00", high = "#FF6600", trans = "reverse") +
#   scale_x_continuous(
#     expand = c(0,0),
#     limits = c(0,14)
#   ) +
#   geom_col(width = .2) +
#   theme(
#     axis.text.y = element_text(color=y_col),
#     axis.text.x = element_text(size = 10)
#   )
ms.GO_CC.df <- ms.GO_CC.df %>%
  dplyr::arrange(desc(Count), Description) %>%
  dplyr::mutate(Description = factor(Description, levels = rev(Description))) %>%
  dplyr::mutate(Count = as.integer(Count))
  
ggplot(ms.GO_CC.df, aes(x = Description, y = -log10(p.adjust), fill = Count)) +
  coord_flip() +
  ylab(bquote(-log[10](p.adj))) +
  geom_bar(stat = "identity", color = "#333333", width = 0.67) +
  scale_y_continuous(
    expand = c(0,0),
    limits = c(0,1.5)
  ) +
  scale_fill_continuous(
    low = "brown", 
    high = "orange",
    limits = c(0,12),
    breaks = c(0,4,8,12)
    ) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 10),
        legend.text = element_text(size = 10)
  )

ggsave(
  filename = "GO_CC_barplot_proteomics.eps", 
  device = "eps", 
  path = "../figures/Fig3/", 
  width = 120, height = 90, 
  dpi = 300, 
  units = "mm", 
  family = "Arial"
)