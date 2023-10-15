library(Seurat)
library(SeuratDisk)
library(RColorBrewer)
library(ggpubr)

suppressMessages(
  extrafont::loadfonts(device="postscript")
)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

hpcs.lps <- LoadH5Seurat("../data/processed/hpcs_lps_state_marked.h5Seurat", verbose = F)

# Control group factors and plot orders
hpcs.lps$state <- factor(hpcs.lps$state, levels = c("Healthy","Inflammation"))
hpcs.lps$cell_type_brief <- factor(hpcs.lps$cell_type_brief, levels = c("Som","Lac","Cort","Mel","Gonad","Thyro"))
hpcs.lps$cell_type_state <- paste(hpcs.lps$cell_type_brief, hpcs.lps$state, sep = "_")
hpcs.lps$cell_type_state <- factor(
  hpcs.lps$cell_type_state, 
  levels = c("Som_Healthy","Som_Inflammation","Lac_Healthy","Lac_Inflammation","Cort_Healthy","Cort_Inflammation",
             "Mel_Healthy","Mel_Inflammation","Gonad_Healthy","Gonad_Inflammation","Thyro_Healthy","Thyro_Inflammation")
    )

plot_violins <- function(object, feature, path){
  object |>
    FetchData(vars = c(feature, "cell_type_state", "cell_type_brief", "state")) |>
    ggviolin(
      x = "cell_type_state",
      y = feature,
      ylab = feature,
      size = .5,
      color = "state",
      fill = "cell_type_state",
      palette = brewer.pal(n = 12, name = "Paired")[c(2,8)],
      width = 1,
      add = list("jitter"),
      add.params = list(color="#333333", size=.5)
    ) +
    scale_y_continuous(limits = c(0,NA)) +
    scale_fill_manual(values = scales::hue_pal()(13)[c(1,1,2,2,9,9,6,6,10,10,13,13)]) +
    theme_pubclean() +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = 'transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 16),
      axis.text.y = element_text(size = 15, color = "black")
    )
  ggsave(
    filename = paste0(feature,"_violin.eps"), device = "eps", path = path, width = 6, height = 2, dpi = 300, family = "Arial"
  )
}

plot_violins(object = hpcs.lps, feature = "Gh", path = "../figures/FigSX/")
plot_violins(object = hpcs.lps, feature = "Prl", path = "../figures/FigSX/")
plot_violins(object = hpcs.lps, feature = "Pomc", path = "../figures/FigSX/")
plot_violins(object = hpcs.lps, feature = "Fshb", path = "../figures/FigSX/")
plot_violins(object = hpcs.lps, feature = "Lhb", path = "../figures/FigSX/")
plot_violins(object = hpcs.lps, feature = "Tshb", path = "../figures/FigSX/")
plot_violins(object = hpcs.lps, feature = "Cga", path = "../figures/FigSX/")

plot_violins(object = hpcs.lps, feature = "Cxcl1", path = "../figures/FigSX/")
plot_violins(object = hpcs.lps, feature = "Cxcl10", path = "../figures/FigSX/")
plot_violins(object = hpcs.lps, feature = "Ccl2", path = "../figures/FigSX/")
plot_violins(object = hpcs.lps, feature = "Ccl5", path = "../figures/FigSX/")
plot_violins(object = hpcs.lps, feature = "Ccl8", path = "../figures/FigSX/")

plot_violins(object = hpcs.lps, feature = "Cartpt", path = "../figures/FigSX/")
plot_violins(object = hpcs.lps, feature = "Fgg", path = "../figures/FigSX/")
plot_violins(object = hpcs.lps, feature = "Fndc4", path = "../figures/FigSX/")
plot_violins(object = hpcs.lps, feature = "Il1r2", path = "../figures/FigSX/")
plot_violins(object = hpcs.lps, feature = "Insl6", path = "../figures/FigSX/")
plot_violins(object = hpcs.lps, feature = "Isg15", path = "../figures/FigSX/")
plot_violins(object = hpcs.lps, feature = "Lbp", path = "../figures/FigSX/")
plot_violins(object = hpcs.lps, feature = "Lgals9", path = "../figures/FigSX/")
plot_violins(object = hpcs.lps, feature = "Lox", path = "../figures/FigSX/")
plot_violins(object = hpcs.lps, feature = "Nptx2", path = "../figures/FigSX/")
plot_violins(object = hpcs.lps, feature = "Nrtn", path = "../figures/FigSX/")
plot_violins(object = hpcs.lps, feature = "Serpina3n", path = "../figures/FigSX/")
plot_violins(object = hpcs.lps, feature = "Vgf", path = "../figures/FigSX/")


x = hpcs.lps |> subset(subset = cell_type_brief == "Lac" & state == "Healthy") |> FetchData(vars = c("Ccl8")) |> pull("Ccl8")
y = hpcs.lps |> subset(subset = cell_type_brief == "Lac" & state == "Inflammation") |> FetchData(vars = c("Ccl8")) |> pull("Ccl8")
wilcox.test(x,y,alternative = "less")

cell_types = c("Som","Lac","Cort","Mel","Gonad","Thyro")
sapply(X = cell_types, FUN = function(x){
  
})

# hpcs.lps |>
#   FetchData(vars = c("Gh", "cell_type_state", "cell_type_brief", "state")) |>
#   ggviolin(
#     x = "cell_type_state",
#     y = "Gh",
#     ylab = "Gh",
#     size = .6,
#     color = "state",
#     # fill = "cell_type_state",
#     palette = brewer.pal(n = 12, name = "Paired")[c(2,8)],
#     width = 1,
#     add = list("jitter"),
#     add.params = list(color="#333333", size=.5)
#   ) +
#   scale_y_continuous(limits = c(0,NA)) +
#   scale_fill_manual(values = scales::hue_pal()(13)[c(1,1,2,2,9,9,6,6,10,10,13,13)]) +
#   theme_pubclean() +
#   theme(
#     legend.position = "none",
#     axis.title.x = element_blank(),
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.title.y = element_text(size = 16),
#     axis.text.y = element_text(size = 15)
#   )

# # Somatotropes: Gh
# som <- subset(hpcs.lps, subset = cell_type_brief == "Som")
# som$state <- factor(som$state, levels = c("Healthy","Inflammation"))
# som |>
#   FetchData(vars = c("Gh","state")) |>
#   ggviolin(
#     x = "state",
#     y = "Gh",
#     ylab = "Normalized expression",
#     title = "Gh",
#     size = .6,
#     fill = "state",
#     palette = brewer.pal(n = 12, name = "Paired")[c(2,8)],
#     width = .7,
#     add = list("jitter"),
#     add.params = list(color="#333333", size=.5)
#   ) +
#   theme_pubr() +
#   theme(
#     legend.position = "none",
#     plot.title = element_text(size = 16),
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(size = 12),
#     axis.title.y = element_text(size = 16)
#   )
# ggsave(filename = "som_Gh_vlnplot.eps", device = "eps", path = "../figures/FigSX/", dpi = 300, width = 5, height = 4)
# 
# 
# # Corticotropes: Pomc
# cort <- subset(hpcs.lps, subset = cell_type_brief == "Cort")
# cort$state <- factor(cort$state, levels = c("Healthy","Inflammation"))
# cort |>
#   FetchData(vars = c("Pomc","state")) |>
#   ggviolin(
#     x = "state",
#     y = "Pomc",
#     ylab = "Normalized expression",
#     title = "Pomc (Corticotropes)",
#     size = .6,
#     fill = "state",
#     palette = brewer.pal(n = 12, name = "Paired")[c(2,8)],
#     width = .7,
#     add = list("jitter"),
#     add.params = list(color="#333333", size=.5)
#   ) +
#   theme_pubr() +
#   theme(
#     legend.position = "none",
#     plot.title = element_text(size = 16),
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(size = 12),
#     axis.title.y = element_text(size = 16)
#   )
# ggsave(filename = "cort_Pomc_vlnplot.eps", device = "eps", path = "../figures/FigSX/", dpi = 300, width = 5, height = 4)
# 
# # Lactotropes: Prl
# lac <- subset(hpcs.lps, subset = cell_type_brief == "Lac")
# lac$state <- factor(lac$state, levels = c("Healthy","Inflammation"))
# lac |>
#   FetchData(vars = c("Prl","state")) |>
#   ggviolin(
#     x = "state",
#     y = "Prl",
#     ylab = "Normalized expression",
#     title = "Prl",
#     size = .6,
#     fill = "state",
#     palette = brewer.pal(n = 12, name = "Paired")[c(2,8)],
#     width = .7,
#     add = list("jitter"),
#     add.params = list(color="#333333", size=.5)
#   ) +
#   theme_pubr() +
#   theme(
#     legend.position = "none",
#     plot.title = element_text(size = 16),
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(size = 12),
#     axis.title.y = element_text(size = 16)
#   )
# ggsave(filename = "lac_Prl_vlnplot.eps", device = "eps", path = "../figures/FigSX/", dpi = 300, width = 5, height = 4)
# 
# # Gonadotropes: Fshb, Lhb, Cga
# gonad <- subset(hpcs.lps, subset = cell_type_brief == "Gonad")
# gonad$state <- factor(gonad$state, levels = c("Healthy","Inflammation"))
# gonad |>
#   FetchData(vars = c("Fshb","state")) |>
#   ggviolin(
#     x = "state",
#     y = "Fshb",
#     ylab = "Normalized expression",
#     title = "Fshb",
#     size = .6,
#     fill = "state",
#     palette = brewer.pal(n = 12, name = "Paired")[c(2,8)],
#     width = .7,
#     add = list("jitter"),
#     add.params = list(color="#333333", size=.5)
#   ) +
#   theme_pubr() +
#   theme(
#     legend.position = "none",
#     plot.title = element_text(size = 16),
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(size = 12),
#     axis.title.y = element_text(size = 16)
#   )
# ggsave(filename = "gonad_Fshb_vlnplot.eps", device = "eps", path = "../figures/FigSX/", dpi = 300, width = 5, height = 4)
# 
# gonad |>
#   FetchData(vars = c("Lhb","state")) |>
#   ggviolin(
#     x = "state",
#     y = "Lhb",
#     ylab = "Normalized expression",
#     title = "Lhb",
#     size = .6,
#     fill = "state",
#     palette = brewer.pal(n = 12, name = "Paired")[c(2,8)],
#     width = .7,
#     add = list("jitter"),
#     add.params = list(color="#333333", size=.5)
#   ) +
#   theme_pubr() +
#   theme(
#     legend.position = "none",
#     plot.title = element_text(size = 16),
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(size = 12),
#     axis.title.y = element_text(size = 16)
#   )
# ggsave(filename = "gonad_Lhb_vlnplot.eps", device = "eps", path = "../figures/FigSX/", dpi = 300, width = 5, height = 4)
# 
# gonad |>
#   FetchData(vars = c("Cga","state")) |>
#   ggviolin(
#     x = "state",
#     y = "Cga",
#     ylab = "Normalized expression",
#     title = "Cga (Gonadotropes)",
#     size = .6,
#     fill = "state",
#     palette = brewer.pal(n = 12, name = "Paired")[c(2,8)],
#     width = .7,
#     add = list("jitter"),
#     add.params = list(color="#333333", size=.5)
#   ) +
#   theme_pubr() +
#   theme(
#     legend.position = "none",
#     plot.title = element_text(size = 16),
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(size = 12),
#     axis.title.y = element_text(size = 16)
#   )
# ggsave(filename = "gonad_Cga_vlnplot.eps", device = "eps", path = "../figures/FigSX/", dpi = 300, width = 5, height = 4)
# 
# # Thyrotropes: Tshb, Cga
# thyro <- subset(hpcs.lps, subset = cell_type_brief == "Thyro")
# thyro$state <- factor(thyro$state, levels = c("Healthy","Inflammation"))
# thyro |>
#   FetchData(vars = c("Tshb","state")) |>
#   ggviolin(
#     x = "state",
#     y = "Tshb",
#     ylab = "Normalized expression",
#     title = "Tshb",
#     size = .6,
#     fill = "state",
#     palette = brewer.pal(n = 12, name = "Paired")[c(2,8)],
#     width = .7,
#     add = list("jitter"),
#     add.params = list(color="#333333", size=.5)
#   ) +
#   theme_pubr() +
#   theme(
#     legend.position = "none",
#     plot.title = element_text(size = 16),
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(size = 12),
#     axis.title.y = element_text(size = 16)
#   )
# ggsave(filename = "thyro_Tshb_vlnplot.eps", device = "eps", path = "../figures/FigSX/", dpi = 300, width = 5, height = 4)
# 
# thyro |>
#   FetchData(vars = c("Cga","state")) |>
#   ggviolin(
#     x = "state",
#     y = "Cga",
#     ylab = "Normalized expression",
#     title = "Cga (Thyrotropes)",
#     size = .6,
#     fill = "state",
#     palette = brewer.pal(n = 12, name = "Paired")[c(2,8)],
#     width = .7,
#     add = list("jitter"),
#     add.params = list(color="#333333", size=.5)
#   ) +
#   theme_pubr() +
#   theme(
#     legend.position = "none",
#     plot.title = element_text(size = 16),
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(size = 12),
#     axis.title.y = element_text(size = 16)
#   )
# ggsave(filename = "thyro_Cga_vlnplot.eps", device = "eps", path = "../figures/FigSX/", dpi = 300, width = 5, height = 4)
# 
# # Melanotropes: Pomc
# mel <- subset(hpcs.lps, subset = cell_type_brief == "Mel")
# mel$state <- factor(mel$state, levels = c("Healthy","Inflammation"))
# mel |>
#   FetchData(vars = c("Pomc","state")) |>
#   ggviolin(
#     x = "state",
#     y = "Pomc",
#     ylab = "Normalized expression",
#     title = "Pomc (Melanotropes)",
#     size = .6,
#     fill = "state",
#     palette = brewer.pal(n = 12, name = "Paired")[c(2,8)],
#     width = .7,
#     add = list("jitter"),
#     add.params = list(color="#333333", size=.5)
#   ) +
#   theme_pubr() +
#   theme(
#     legend.position = "none",
#     plot.title = element_text(size = 16),
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(size = 12),
#     axis.title.y = element_text(size = 16)
#   )
# ggsave(filename = "mel_Pomc_vlnplot.eps", device = "eps", path = "../figures/FigSX/", dpi = 300, width = 5, height = 4)
# 
# obj.list = list(som=som, cort=cort, lac=lac, gonad=gonad, thyro=thyro, mel=mel)
# for (obj in obj.list){
#   obj
# }
