library(Seurat)
library(SeuratDisk)
library(tidyverse)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
suppressMessages(
  extrafont::loadfonts(device="postscript")
)

source("utilities.R")


hpcs.lps <- LoadH5Seurat("../data/processed/hpcs_lps_state_marked.h5Seurat")

tab_xiong <- readxl::read_excel("../misc/Xiong_et_al_mouse_secretome_receptor_list_curated.xlsx")

GO_CC <- read.csv("../outs/GO_CC.csv", row.names = 1)
# go_features <- GO_CC %>% dplyr::filter(Description == "extracellular space") %>% pull(geneID) %>% str_split("/") %>% unlist()
go_features <- GO_CC %>% dplyr::filter(Description == "extracellular region") %>% pull(geneID) %>% str_split("/") %>% unlist()
# go_features <- go_features[go_features %in% c("H2-D1","")]  ## Some of these candidates doesn't conform to prior biological knowledge

test_method = "MAST"
de.state.mks <- read.csv(file = paste0("../outs/de_state_markers_",test_method,".csv"))

features <- union(go_features, tab_xiong %>% dplyr::filter(Type == "secreted") %>% pull("Symbol"))
features <- features[features %in% de.state.mks$gene]
features <- sort(c(features, "Lbp"))

hpcs.lps$cell_type_brief <- factor(hpcs.lps$cell_type_brief, levels = c("Som","Lac","Cort","Mel","Gonad","Thyro"))
hpcs.lps$state <- factor(hpcs.lps$state, levels = c("Healthy","Inflammation"))

hpcs.lps <- ScaleData(hpcs.lps, features = c(VariableFeatures(hpcs.lps), features))
hpcs.lps$cell_type_brief <- factor(hpcs.lps$cell_type_brief, levels = c("Som","Lac","Cort","Mel","Gonad","Thyro"))
hpcs.lps$state <- factor(hpcs.lps$state, levels = c("Healthy","Inflammation"))
# hpcs.lps <- PercentExpr(hpcs.lps)

hpcs.lps.healthy <- subset(hpcs.lps, subset = state == "Healthy")
hpcs.lps.inflammation <- subset(hpcs.lps, subset = state == "Inflammation")

### Healthy
hpcs.lps.healthy <- CalculateMarkerSpecificity(hpcs.lps.healthy, group.by = "cell_type_brief")
df.spec.healthy <- hpcs.lps.healthy[["RNA"]][[]] %>% dplyr::select(starts_with("spec."))
df.spec.healthy <- df.spec.healthy %>% rownames_to_column("gene") %>%
  dplyr::filter(gene %in% features)
df.spec.healthy <- df.spec.healthy %>% 
  pivot_longer(
    cols = starts_with("spec."),
    names_to = "cell_type",
    values_to = "specificity"
  ) %>%
  mutate(cell_type = sub("spec.","",cell_type))

df.avgexp.healthy <- AverageExpression(object = hpcs.lps.healthy, assays = "RNA", features = features, group.by = "cell_type_brief")
df.avgexp.healthy <- df.avgexp.healthy[[1]] %>% as.data.frame() %>%
  # rename_with(.fn = ~ paste0("avgexp.", .x)) %>%
  rownames_to_column("gene")

df.avgexp.healthy <- df.avgexp.healthy %>%
  pivot_longer(
    cols = c(Som,Lac,Cort,Mel,Gonad,Thyro),
    names_to = "cell_type",
    values_to = "avg.expr"
  ) 
# %>%  mutate(cell_type = sub("avgexp.","",cell_type))

hpcs.lps.heal.list <- SplitObject(object = hpcs.lps.healthy, split.by = "cell_type_brief")
hpcs.lps.heal.list <- lapply(hpcs.lps.heal.list, PercentExpr)

df.pct.expr.healthy <- data.frame(row.names = features)
for (i in seq_along(hpcs.lps.heal.list)){
  cell_type = names(hpcs.lps.heal.list[i])
  df.pct.expr.healthy[[cell_type]] <- hpcs.lps.heal.list[[i]][["RNA"]][[]][features,"pct.expr"]
}

df.pct.expr.healthy <- df.pct.expr.healthy %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = c(Som,Lac,Cort,Mel,Gonad,Thyro),
    names_to = "cell_type",
    values_to = "pct.expr"
  )

df.healthy <- df.spec.healthy %>% left_join(df.avgexp.healthy, by = c("gene","cell_type")) %>% left_join(df.pct.expr.healthy, by = c("cell_type", "gene"))
df.healthy$state <- "Healthy"

### Inflammation
hpcs.lps.inflammation <- CalculateMarkerSpecificity(hpcs.lps.inflammation, group.by = "cell_type_brief")
df.spec.inflammation <- hpcs.lps.inflammation[["RNA"]][[]] %>% dplyr::select(starts_with("spec."))
df.spec.inflammation <- df.spec.inflammation %>% rownames_to_column("gene") %>%
  dplyr::filter(gene %in% features)
df.spec.inflammation <- df.spec.inflammation %>% 
  pivot_longer(
    cols = starts_with("spec."),
    names_to = "cell_type",
    values_to = "specificity"
  ) %>%
  mutate(cell_type = sub("spec.","",cell_type))

df.avgexp.inflammation <- AverageExpression(object = hpcs.lps.inflammation, assays = "RNA", features = features, group.by = "cell_type_brief")
df.avgexp.inflammation <- df.avgexp.inflammation[[1]] %>% as.data.frame() %>%
  # rename_with(.fn = ~ paste0("avgexp.", .x)) %>%
  rownames_to_column("gene")

df.avgexp.inflammation <- df.avgexp.inflammation %>%
  pivot_longer(
    cols = c(Som,Lac,Cort,Mel,Gonad,Thyro),
    names_to = "cell_type",
    values_to = "avg.expr"
  ) 
# %>% mutate(cell_type = sub("avgexp.","",cell_type))

hpcs.lps.infl.list <- SplitObject(object = hpcs.lps.inflammation, split.by = "cell_type_brief")
hpcs.lps.infl.list <- lapply(hpcs.lps.infl.list, PercentExpr)

df.pct.expr.inflammation <- data.frame(row.names = features)
for (i in seq_along(hpcs.lps.infl.list)){
  cell_type = names(hpcs.lps.infl.list[i])
  df.pct.expr.inflammation[[cell_type]] <- hpcs.lps.infl.list[[i]][["RNA"]][[]][features,"pct.expr"]
}

df.pct.expr.inflammation <- df.pct.expr.inflammation %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = c(Som,Lac,Cort,Mel,Gonad,Thyro),
    names_to = "cell_type",
    values_to = "pct.expr"
  )

df.inflammation <- df.spec.inflammation %>% left_join(df.avgexp.inflammation, by = c("gene","cell_type")) %>% left_join(df.pct.expr.inflammation, by = c("cell_type", "gene"))
df.inflammation$state <- "Inflamamtion"

df.healthy$cell_type <- factor(df.healthy$cell_type, levels = rev(c("Som","Lac","Cort","Mel","Gonad","Thyro")))
df.inflammation$cell_type <- factor(df.inflammation$cell_type, levels = rev(c("Som","Lac","Cort","Mel","Gonad","Thyro")))

features_qpcr = c("Cartpt","Nptx2","Lgals9","Nrtn","Serpina3n","Lox","Insl6","Fndc4","Lbp","Fgg","Isg15","Vgf","Il1r2")

d1 <- df.healthy %>% 
  dplyr::filter(gene %in% features_qpcr) %>%
  ggplot(mapping = aes(x = gene, y = cell_type)) +
  geom_point(mapping = aes(size = pct.expr, color = specificity)) +
  scale_size(limits = c(0,1), breaks = c(.2,.4,.6)) +
  scale_color_distiller(palette = "YlGnBu", direction = 1) +
  theme_test()

d2 <- df.inflammation %>% 
  dplyr::filter(gene %in% features_qpcr) %>%
  ggplot(mapping = aes(x = gene, y = cell_type)) +
  geom_point(mapping = aes(size = pct.expr, color = specificity)) +
  labs(y = "Pituitary Cell Type") +
  scale_size(limits = c(0,1), breaks = c(.2,.4,.6)) +
  scale_color_distiller(palette = "YlOrRd", direction = 1) +
  theme_test() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  )
ggsave(
  filename = "dotplot_infl_spec_pctexpr.eps",
  plot = d2,
  device = "eps",
  path = "../figures/Fig3/",
  width = 12, height = 3.6,
  dpi = 300,
  family = "Arial"
)

d2_clean <- df.inflammation %>% 
  dplyr::filter(gene %in% features_qpcr) %>%
  ggplot(mapping = aes(x = gene, y = cell_type)) +
  geom_point(mapping = aes(size = pct.expr, color = specificity)) +
  scale_size(limits = c(0,1), breaks = c(.2,.4,.6)) +
  scale_color_distiller(palette = "YlOrRd", direction = 1) +
  theme_test() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 11, face = "bold")
  )
ggsave(
  filename = "dotplot_infl_spec_pctexpr_clean.eps",
  plot = d2_clean,
  device = "eps",
  path = "../figures/Fig3/",
  width = 10, height = 2.5,
  dpi = 300,
  family = "Arial"
)