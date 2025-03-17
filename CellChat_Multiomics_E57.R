
#This is the code for E57 Multiomics CellChat analysis. 

load("E57.RData")

# For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
labels <- Idents(IL)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

data.input = Seurat::GetAssayData(IL, slot = "data", assay = "SCT") # normalized data matrix
meta = data.frame(labels = Seurat::Idents(IL), samples = "sample1", row.names = names(Seurat::Idents(IL))) # manually create a dataframe consisting of the cell labels


library(CellChat)
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

cellChat <- createCellChat(object = IL, group.by = "ident", assay = "RNA")


CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)



# use Secreted Signaling
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") 

cellchat@DB <- CellChatDB.use



# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 692

#Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, type = "triMean")


cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)


#cell chat interactions between cell types and interactions
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


mat <- cellchat@net$weight
par(mfrow = c(3,4))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

#outgoing and incoming signaling patterns for all cell types
#based on this, you can do subsequent network analysis
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

#To look at outgoing communication patterns and secreting cells
library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
dev.off()
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
netAnalysis_river(cellchat, pattern = "outgoing")



#Cell types and interactions with other cell types
#Comparing communications on a single object, # (1) show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')

#RestingChondro
netVisual_bubble(cellchat, sources.use = 8, targets.use = c(8,10,13,17), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 8, targets.use = c(1:19), remove.isolate = FALSE)
#perichondro
netVisual_bubble(cellchat, sources.use = 13, targets.use = c(8,10,13,17), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 13, targets.use = c(1:19), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 10, targets.use = c(1:10), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:10), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 2, targets.use = c(1:10), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:10), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:10), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:10), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:10), remove.isolate = FALSE)


cellchat@netP$pathways
#Wnt, TGF-β, BMP, FGF, and Hedgehog 

pathways.show <- c("ncWNT") 
pathways.show <- c("PTN") 
pathways.show <- c("WNT")
pathways.show <- c("BMP") 
pathways.show <- c("SEMA3") 
pathways.show <- c("GDF") 
pathways.show <- c("MK") 
 

pathways.show <- c("FGF") 
pathways.show <- c("PTH")
pathways.show <- c("TGFb")
pathways.show <- c("VEGF")
pathways.show <- c("PTPR")
pathways.show <- c("IGF")
pathways.show <- c("IGFBP") 
pathways.show <- c("PTPR") 
pathways.show <- c("PERIOSTIN")
pathways.show <- c("HH")
pathways.show <- c("SLIT")
pathways.show <- c("GDF")
pathways.show <- c("LIFR")
pathways.show <- c("SLITRK")


# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

#>Do heatmap based on a single object

#Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
netAnalysis_contribution(cellchat, signaling = pathways.show)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 15, height = 4.5, font.size = 10)


