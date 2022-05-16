library(Seurat)
library(SeuratDisk)
library(tidyverse)
source("utilities.R")

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

hpcs.lps <- LoadH5Seurat("../data/processed/hpcs_lps_state_marked.h5Seurat")

tab_xiong <- readxl::read_excel("../misc/Xiong_et_al_mouse_secretome_receptor_list_curated.xlsx")

GO_CC <- read.csv("../outs/GO_CC.csv", row.names = 1)
go_features <- GO_CC %>% dplyr::filter(Description == "extracellular space") %>% pull(geneID) %>% str_split("/") %>% unlist()
# go_features <- go_features[go_features %in% c("H2-D1","")]  ## Some of these candidates doesn't conform to prior biological knowledge

features <- union(go_features, tab_xiong %>% dplyr::filter(Type == "secreted") %>% pull("Symbol"))

hpcs.lps$cell_type_brief <- factor(hpcs.lps$cell_type_brief, levels = c("Som","Lac","Cort","Mel","Gonad","Thyro"))
hpcs.lps$state <- factor(hpcs.lps$state, levels = c("Healthy","Inflammation"))

hpcs.lps <- ScaleData(hpcs.lps, features = c(VariableFeatures(hpcs.lps), features))

hpcs.lps <- PercentExpr(hpcs.lps)

hpcs.lps.healthy <- subset(hpcs.lps, subset = state == "Healthy")
hpcs.lps.inflammation <- subset(hpcs.lps, subset = state == "Inflammation")

### Healthy
hpcs.lps.healthy <- CalculateMarkerSpecificity(hpcs.lps.healthy, idents = "cell_type_brief")
df.avgexp.healthy <- AverageExpression(object = hpcs.lps.healthy, assays = "RNA", features = features, group.by = "cell_type_brief")
df.avgexp.healthy <- df.avgexp.healthy[[1]] %>% as.data.frame() %>% rename_with(.fn = ~ paste0("avgexp.", .x)) %>% rownames_to_column("gene")
df.avgexp.healthy <- df.avgexp.healthy %>%
  pivot_longer(
    cols = starts_with("avgexp."),
    names_to = "cell_type",
    values_to = "avg.expr"
  ) %>%
  mutate(cell_type = sub("avgexp.","",cell_type))

df.spec.healthy <- hpcs.lps.healthy[["RNA"]][[]] %>% dplyr::select(pct.expr, starts_with("spec."))
df.spec.healthy <- df.spec.healthy %>% rownames_to_column("gene") %>%
  dplyr::filter(gene %in% features)
df.spec.healthy <- df.spec.healthy %>% 
  pivot_longer(
    cols = starts_with("spec."),
    names_to = "cell_type",
    values_to = "specificity"
  ) %>%
  mutate(cell_type = sub("spec.","",cell_type))

df.healthy <- left_join(x = df.spec.healthy, y = df.avgexp.healthy, by = c("gene","cell_type"))

### Inflammation
hpcs.lps.inflammation <- CalculateMarkerSpecificity(hpcs.lps.inflammation, idents = "cell_type_brief")
df.avgexp.inflammation <- AverageExpression(object = hpcs.lps.inflammation, assays = "RNA", features = features, group.by = "cell_type_brief")
df.avgexp.inflammation <- df.avgexp.inflammation[[1]] %>% as.data.frame() %>% rename_with(.fn = ~ paste0("avgexp.", .x)) %>% rownames_to_column("gene")
df.avgexp.inflammation <- df.avgexp.inflammation %>%
  pivot_longer(
    cols = starts_with("avgexp."),
    names_to = "cell_type",
    values_to = "avg.expr"
  ) %>%
  mutate(cell_type = sub("avgexp.","",cell_type))

df.spec.inflammation <- hpcs.lps.inflammation[["RNA"]][[]] %>% dplyr::select(pct.expr, starts_with("spec."))
df.spec.inflammation <- df.spec.inflammation %>% rownames_to_column("gene") %>%
  dplyr::filter(gene %in% features)
df.spec.inflammation <- df.spec.inflammation %>% 
  pivot_longer(
    cols = starts_with("spec."),
    names_to = "cell_type",
    values_to = "specificity"
  ) %>%
  mutate(cell_type = sub("spec.","",cell_type))

df.inflammation <- left_join(x = df.spec.inflammation, y = df.avgexp.inflammation, by = c("gene","cell_type"))
