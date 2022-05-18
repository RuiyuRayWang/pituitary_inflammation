library(Seurat)
library(SeuratDisk)
library(CellChat)
library(patchwork)
library(tidyverse)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

suppressMessages(
  extrafont::loadfonts(device="postscript")
)

# Part I: Data input & processing and initialization of CellChat object
## Load data
hpcs.lps <- LoadH5Seurat("../data/processed/hpcs_lps_state_marked.h5Seurat")
hpcs.lps.inflammation <- subset(hpcs.lps, subset = state == "Inflammation")

### Normalized data (e.g., library-size normalization and then log-transformed with a pseudocount of 1) is required as input for CellChat analysis.
data.input = GetAssayData(hpcs.lps.inflammation, slot = "data", assay = "RNA")
meta <- hpcs.lps.inflammation@meta.data

## Create a CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type_brief")

## Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

### Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

### use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

### set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 16) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# # Optional: project gene expression data onto PPI (when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)


# Part II: Inference of cell-cell communication network
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

### Extract the inferred cellular communication network as a data frame
# returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. 
# Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
df.net <- subsetCommunication(cellchat)
# df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))  # inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))  # inferred cell-cell communications mediated by signaling WNT and TGFb.

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)


# Part III: Visualization of cell-cell communication network



# Part IV: Systems analysis of cell-cell communication network

## Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling

### Compute and visualize the network centrality scores
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways


# Part V: Save the CellChat object
saveRDS(cellchat, file = "../data/cellchat/cellchat_hpcslps_inflammation.rds")
