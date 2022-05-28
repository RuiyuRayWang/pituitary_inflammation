library(tidyverse)
library(nVennR)

test_method = "MAST"
csvd.cell.mks <- read.csv(file = paste0("../outs/conserved_cell_markers_",test_method,".csv"))

csvd_V <- plotVenn(
  sets = list(
    Som=csvd.cell.mks %>% dplyr::filter(cell_type == "Som") %>% pull("gene"), 
    Lac=csvd.cell.mks %>% dplyr::filter(cell_type == "Lac") %>% pull("gene"), 
    Cort=csvd.cell.mks %>% dplyr::filter(cell_type == "Cort") %>% pull("gene"),
    Mel=csvd.cell.mks %>% dplyr::filter(cell_type == "Mel") %>% pull("gene"),
    Gonad=csvd.cell.mks %>% dplyr::filter(cell_type == "Gonad") %>% pull("gene"),
    Thyro=csvd.cell.mks %>% dplyr::filter(cell_type == "Thyro") %>% pull("gene")
    ),
  nCycles = 10000,
  sNames = list("Som","Lac","Cort","Mel","Gonad","Thyro"),
  showPlot = F
  )

showSVG(
  csvd_V, 
  opacity=0.3, 
  fontScale = 2,
  setColors = c("#F8766D", "#E18A00", "#00ACFC", "#00BE70", "#8B93FF", "#FF65AC"),
  outFile = "../figures/Fig2/nVenn_csvd.svg"
)
showSVG(
  csvd_V, 
  showLegend = F,
  labelRegions = F,
  opacity=0.3, 
  fontScale = 2,
  setColors = c("#F8766D", "#E18A00", "#00ACFC", "#00BE70", "#8B93FF", "#FF65AC"),
  outFile = "../figures/Fig2/nVenn_csvd_clean.svg"
  )

de.state.mks <- read.csv(file = paste0("../outs/de_state_markers_",test_method,".csv"))

de_V <- plotVenn(
  sets = list(
    Som=de.state.mks %>% dplyr::filter(cell_type == "Som") %>% pull("gene"), 
    Lac=de.state.mks %>% dplyr::filter(cell_type == "Lac") %>% pull("gene"), 
    Cort=de.state.mks %>% dplyr::filter(cell_type == "Cort") %>% pull("gene"),
    Mel=de.state.mks %>% dplyr::filter(cell_type == "Mel") %>% pull("gene"),
    Gonad=de.state.mks %>% dplyr::filter(cell_type == "Gonad") %>% pull("gene"),
    Thyro=de.state.mks %>% dplyr::filter(cell_type == "Thyro") %>% pull("gene")
  ),
  nCycles = 17000,
  sNames = list("Som","Lac","Cort","Mel","Gonad","Thyro"),
  showPlot = F
)
showSVG(
  de_V, 
  opacity=0.3, 
  fontScale = 2,
  setColors = c("#F8766D", "#E18A00", "#00ACFC", "#00BE70", "#8B93FF", "#FF65AC"),
  outFile = "../figures/Fig2/nVenn_destate.svg"
)
showSVG(
  de_V, 
  showLegend = F,
  labelRegions = F,
  opacity=0.3, 
  fontScale = 2,
  setColors = c("#F8766D", "#E18A00", "#00ACFC", "#00BE70", "#8B93FF", "#FF65AC"),
  outFile = "../figures/Fig2/nVenn_destate_clean.svg"
  )