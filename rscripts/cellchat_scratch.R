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

group.cellType.saline <- c("Spleen","Pituitary","Pituitary","Pituitary","Spleen","Pituitary","Spleen","Spleen","Spleen","Spleen","Pituitary","Spleen","Pituitary")
names(group.cellType.saline) <- levels(cellchat.saline@idents)
group.cellType.saline
pit.saline.group.idx <- which(group.cellType.saline == "Pituitary")
spleen.saline.group.idx <- which(group.cellType.saline == "Spleen")

group.cellType.lps <- c("Spleen","Pituitary","Pituitary","Pituitary","Spleen","Pituitary","Spleen","Spleen","Spleen","Spleen","Pituitary","Spleen","Pituitary")
names(group.cellType.lps) <- levels(cellchat.lps@idents)
group.cellType.saline
pit.saline.group.idx <- which(group.cellType.saline == "Pituitary")
spleen.saline.group.idx <- which(group.cellType.saline == "Spleen")

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


pathways.show.all <- cellchat@netP$pathways

cellchat@idents <- factor(x = cellchat@idents, levels = c(
  "B cells", "T cells", "NK cells", "pDCs", 
  "Somatotropes", "Lactotropes", "Corticotropes", "Melanotropes", "Gonadotropes", "Thyrotropes", 
  "Monocytes", "Macrophages", "Neutrophils"
))
vertex.receiver = c(5,6,7,8,9,0)
par(mfrow = c(1,1), xpd=TRUE)
for (i in 1:length(pathways.show.all)){
  netVisual_aggregate(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "circle")
  netVisual_aggregate(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "circle", cell.order = c(
    "B cells", "T cells", "NK cells", "pDCs", 
    "Somatotropes", "Lactotropes", "Corticotropes", "Melanotropes", "Gonadotropes", "Thyrotropes", 
    "Monocytes", "Macrophages", "Neutrophils"
  ))
  netVisual_aggregate(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
}






