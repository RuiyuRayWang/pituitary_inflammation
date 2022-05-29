library(CellChat)
library(patchwork)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

suppressMessages(
  extrafont::loadfonts(device="postscript")
)

cellchat.saline <- readRDS("../data/cellchat/cellchat_pituitary_spleen_healthy.rds")
cellchat.lps <- readRDS("../data/cellchat/cellchat_pituitary_spleen_inflammation.rds")
df.net.saline <- subsetCommunication(cellchat.saline)
df.net.lps <- subsetCommunication(cellchat.lps)

group.cellType <- c("Spleen","Pituitary","Pituitary","Pituitary","Spleen","Pituitary","Spleen","Spleen","Spleen","Spleen","Pituitary","Spleen","Pituitary")
names(group.cellType) <- levels(cellchat.saline@idents)
# names(group.cellType) <- c("B","Cort","Gonad","Lac","Macro","Mel","Mono","Neut","NK","pDC","Som","T","Thyro")
group.cellType
pit.group.idx <- which(group.cellType == "Pituitary")
spleen.group.idx <- which(group.cellType == "Spleen")

groupSize <- as.numeric(table(cellchat.saline@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat.saline@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Number of interactions")
netVisual_circle(cellchat.saline@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

groupSize <- as.numeric(table(cellchat.lps@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat.lps@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat.lps@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Automatically save the plots of the all inferred network for quick exploration
dir.create("../outs/cellchat/Saline", recursive = T)
pathways.show.all.saline <- cellchat.saline@netP$pathways
for (i in 1:length(pathways.show.all.saline)){
  try(
    {
      if (any(match(df.net.saline %>% filter(pathway_name == pathways.show.all.saline[i]) %>% pull(source),levels(cellchat.saline@idents)) %in% pit.group.idx) &
          any(match(df.net.saline %>% filter(pathway_name == pathways.show.all.saline[i]) %>% pull(target),levels(cellchat.saline@idents)) %in% spleen.group.idx))
      {
        # Chord diagram
        svglite::svglite(filename = file.path("../outs/cellchat/Saline", paste(pathways.show.all.saline[i],"chort","aggregate.svg", sep = "_")), 
                         width = 6, height = 6)
        netVisual_chord_cell(cellchat.saline, signaling = pathways.show.all.saline[i], group = group.cellType,
                             sources.use = pit.group.idx, targets.use = spleen.group.idx, scale = TRUE)
        dev.off()
      } else {
        print(paste0("Skipping ",pathways.show.all.saline[i]," since no Pituitary-Spleen communication could be found for this pathway."))
      }
    }
  )
}

dir.create("../outs/cellchat/LPS", recursive = T)
pathways.show.all.lps <- cellchat.lps@netP$pathways
for (i in 1:length(pathways.show.all.lps)){
  try(
    {
      if (any(match(df.net.lps %>% filter(pathway_name == pathways.show.all.lps[i]) %>% pull(source),levels(cellchat.lps@idents)) %in% pit.group.idx) &
          any(match(df.net.lps %>% filter(pathway_name == pathways.show.all.lps[i]) %>% pull(target),levels(cellchat.lps@idents)) %in% spleen.group.idx))
      {
        # Chord diagram
        svglite::svglite(filename = file.path("../outs/cellchat/LPS", paste(pathways.show.all.lps[i],"chort","aggregate.svg", sep = "_")), 
                         width = 6, height = 6)
        netVisual_chord_cell(cellchat.lps, signaling = pathways.show.all.lps[i], group = group.cellType, 
                             sources.use = pit.group.idx, targets.use = spleen.group.idx, scale = TRUE)
        dev.off()
      } else {
        print(paste0("Skipping ",pathways.show.all.lps[i]," since no Pituitary-Spleen communication could be found for this pathway."))
      }
    }
  )
}

## Comparative analysis
object.list <- list(Saline = cellchat.saline, LPS = cellchat.lps)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


### Part 1: Whether the cell-cell communication is enhanced or not
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1300)) +
  theme(text = element_text(size = 16))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 100)) +
  theme(text = element_text(size = 16))
gg1 + gg2

### Part 2: The interaction between which cell types is significantly changed
group.cellType <- factor(group.cellType, levels = c("Pituitary","Spleen"))
object.list <- lapply(object.list, function(x){mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents", "counts", "count.merged"))
svglite::svglite(file.path("../outs/cellchat", "interaction_aggregate.svg"), width = 9, height = 4)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1,5,7,8,9,10,12),  comparison = c(1, 2), angle.x = 45)
# ggsave(filename = "comp_cort_pituitary_bubble.eps", device = "eps", plot = last_plot(), width = 6, height = 5, dpi = 300, 
#        path = file.path(out_dir,"comparison"))

# gg1 <- netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1,5,7,8,9,10,12),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LPS", angle.x = 45, remove.isolate = T)
# #> Comparing communications on a merged object
# gg2 <- netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1,5,7,8,9,10,12),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LPS", angle.x = 45, remove.isolate = T)
# #> Comparing communications on a merged object
# gg1 + gg2
# ggsave(filename = "comp_cort_pituitary_bubble_paired.eps", device = "eps", plot = last_plot(), width = 9, height = 5, dpi = 300, 
#        path = file.path(out_dir,"comparison"))

