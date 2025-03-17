#This is the code used for E45_spatial CellChat analysis.

data_dir <- 'E45_01' #ventral section
spatial <- Load10X_Spatial(data_dir, filename = "filtered_feature_bc_matrix.h5", 
                           assay = "Spatial", #specify name of the initial assay
                           slice = "slice1",
                           filter.matrix = TRUE,
                           to.upper = FALSE,
                           image = NULL)
plot1 <- VlnPlot(spatial, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()

plot2 <- SpatialFeaturePlot(spatial, features = "nCount_Spatial",pt.size = 4)
wrap_plots(plot1, plot2)

#SCTransforming the dataset
spatial <- SCTransform(spatial, assay = "Spatial", verbose = FALSE)

#We can then proceed to run dimensionality reduction and clustering on the RNA expression data, using the same workflow as we use for scRNA-seq analysis.
spatial <- RunPCA(spatial, assay = "SCT", verbose = FALSE)
spatial <- FindNeighbors(spatial, reduction = "pca", dims = 1:30)
spatial <- FindClusters(spatial, verbose = FALSE)
spatial <- RunUMAP(spatial, reduction = "pca", dims = 1:30)

Idents(spatial)

spatial <- RenameIdents(spatial, '0' = 'Mes1', '1' = 'Endothelial1','2' = 'Mes2','3' = 'MyoD')
spatial <- RenameIdents(spatial, '4' = 'RestingChondro', '5' ='Chondro1', '6' = 'Fibro1')
spatial <- RenameIdents(spatial, '7' = 'Chondro2', '8' ='ChondroProg', '9' = 'Epithelial')
spatial <- RenameIdents(spatial, '10' = 'Muscles1', '11' ='Mes3', '12' = 'Epithelial2')

spatial$celltype <- Idents(spatial)
Seurat::SpatialDimPlot(spatial)
save(spatial, file = "visium_45_01.RData")

load("visium_45_01.RData")

rna_markers1 <- FindMarkers(spatial, ident.1 = 12, assay = "SCT")
head(rna_markers1)
head(rna_markers1, n = 50)



# show the image and annotated spots
color.use <- scPalette(nlevels(spatial)); names(color.use) <- levels(spatial)


spatial <- RenameIdents(spatial, '0' = 'Mes1', '1' = 'Endothelial1','2' = 'Mes2','3' = 'Mes2')
spatial <- RenameIdents(spatial, '4' = 'RestingChondro', '5' ='Chondro1', '6' = 'Fibro1')
spatial <- RenameIdents(spatial, '7' = 'Chondro2', '8' ='ChondroProg', '9' = 'Epithelial')
spatial <- RenameIdents(spatial, '10' = 'Muscles1', '11' ='Mes3', '12' = 'Epithelial2')

cluster_colors <- c("Mes1" = "#a45c40", "Endothelial1" = "#f0ccb0", "Mes2" = "#523a28", "Mes2"="#523a28",
                    "RestingChondro" = "#145da0", "Chondro1" = "#2e8bc0","Fibro1" = "#59981a","Chondro2" = "#0c2d48","ChondroProg" = "#b1d4e0","Epithelial" = "#a45c40",
                    "Muscles1" = "#ef7c8e","Mes3" = "#d48c70","Epithelial2" = "#a49393")


Seurat::SpatialDimPlot(spatial, label = T, label.size = 3, cols = cluster_colors, pt.size.factor = 3)



# Prepare input data for CelChat analysis
data.input = Seurat::GetAssayData(spatial, slot = "data", assay = "SCT") # normalized data matrix

meta = data.frame(labels = Seurat::Idents(spatial), samples = "sample1", row.names = names(Seurat::Idents(spatial))) # manually create a dataframe consisting of the cell labels
meta$samples <- factor(meta$samples)
unique(meta$labels) # check the cell labels
unique(meta$samples)  # check the sample labels


# load spatial imaging information
# Spatial locations of spots from full (NOT high/low) resolution images are required
spatial.locs = GetTissueCoordinates(spatial, scale = NULL, cols = c("imagerow", "imagecol")) 
# Scale factors and spot diameters of the full resolution images

scale.factors = jsonlite::fromJSON(txt = file.path("/Users/gayanisenevirathne/Desktop/Desktop/Single_cell_Data_Analysis/Ilium_Multiomics_Analysis/Multiomics_E57_IL/outs_57_02_new/spatial", 'scalefactors_json.json'))
scale.factors = list(spot.diameter = 27, spot = scale.factors$spot_diameter_fullres, # these two information are required
                     fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef # these three information are not required
)


coordinates <- as.matrix(spatial.locs[, c("x", "y")])
n_cells <- ncol(data.input)
nrow(coordinates)

n_cells <- nrow(spatial.locs)
spatial.factors <- data.frame(ratio = rep(1, n_cells), tol = rep(0.1, n_cells))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = coordinates, spatial.factors = spatial.factors)


CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)


# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use



# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4)# do parallel

#> Warning: Strategy 'multiprocess' is deprecated in future (>= 1.20.0). Instead,
#> explicitly specify either 'multisession' or 'multicore'. In the current R
#> session, 'multiprocess' equals 'multisession'.
#> Warning in supportsMulticoreAndRStudio(...): [ONE-TIME WARNING] Forked
#> processing ('multicore') is not supported when running R from RStudio
#> because it is considered unstable. For more details, how to control forked
#> processing or not, and how to silence this warning in future R sessions, see ?
#> parallelly::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.mouse)



cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                              distance.use = TRUE, contact.knn.k = 2, scale.distance = 2)

#> truncatedMean is used for calculating the average gene expression per cell group. 
#> [1] ">>> Run CellChat on spatial imaging data using distances as constraints <<< [2022-11-12 07:49:23]"
#> The suggested minimum value of scaled distances is in [1,2], and the calculated value here is  1.30553 
#> [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2022-11-12 08:10:42]"

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)


cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)


groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


library(NMF)
library(ggalluvial)

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

#ChondroProg
netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:13), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 12, targets.use = c(1:13), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 8, targets.use = c(1:13), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 9, targets.use = c(1:13), remove.isolate = FALSE, )
netVisual_bubble(cellchat, sources.use = 6, targets.use = c(1:13), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 13, targets.use = c(1:13), remove.isolate = FALSE)

#to save the graphs as a .pdf use this command
p <- netVisual_bubble(cellchat, sources.use = 9, targets.use = c(1:13), remove.isolate = FALSE, )
p + theme(plot.margin = margin(1, 1, 1, 1, "cm"))  # Adjust margins
ggsave("bubble_plot.pdf", plot = p, width = 10, height = 16, dpi = 300)




cellchat@netP$pathways

pathways.show <- c("MK")

pathways.show <- c("BMP") 
pathways.show <- c("SEMA3") 
pathways.show <- c("GDF") 
pathways.show <- c("PERIOSTIN") 
pathways.show <- c("ncWNT") 
pathways.show <- c("FGF") 
pathways.show <- c("PTH")
pathways.show <- c("PTPR")
pathways.show <- c("EGF")
pathways.show <- c("TGFb")
pathways.show <- c("MK")


pathways.show <- c("PTN")
pathways.show <- c("VEGF")
pathways.show <- c("HH")
pathways.show <- c("ACTIVIN")
pathways.show <- c("WNT")
pathways.show <- c("IGF")
pathways.show <- c("CHEMERIN")
pathways.show <- c("PARs")
pathways.show <- c("CSF")



# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")


# Spatial plot
par(mfrow=c(1,1))

#spatial plot with cluster colors

netVisual_aggregate(cellchat, 
                    signaling = pathways.show, 
                    layout = "spatial", 
                    edge.width.max = 2, 
                    vertex.size.max = 1, 
                    alpha.image = 0.3, 
                    vertex.label.cex = 3.5, 
                    point.size = 4)     # Specify which metadata column to color by


# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 15, height = 4.5, font.size = 10)



# USER can visualize this information on the spatial imaging, e.g., bigger circle indicates larger incoming signaling
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", point.size = 2, edge.width.max = 2, alpha.image = 0.4, vertex.weight = "incoming", vertex.size.max = 3, vertex.label.cex = 3.5)


# Take an input of a few genes
spatialFeaturePlot(cellchat, features = c("PTHLH","PTH1R"), point.size = 1.0, color.heatmap = "Reds", direction = 1)
spatialFeaturePlot(cellchat, features = c("PTN","SDC4"), point.size = 1.0, color.heatmap = "Reds", direction = 1)
spatialFeaturePlot(cellchat, features = c("PTN","SDC2"), point.size = 1.0, color.heatmap = "Reds", direction = 1)
spatialFeaturePlot(cellchat, features = c("MDK","SDC4"), point.size = 1.0, color.heatmap = "Reds", direction = 1)
spatialFeaturePlot(cellchat, features = c("PTN","NCL"), point.size = 2.0, color.heatmap = "Reds", direction = 1)
spatialFeaturePlot(cellchat, features = c("WNT5A","FZD9"), point.size = 1.0, color.heatmap = "Reds", direction = 1)


# Take an input of a ligand-receptor pair
spatialFeaturePlot(cellchat, pairLR.use = "PTHLH_PTH1R", point.size = 1.0, do.binary = FALSE, cutoff = 1.0, enriched.only = F, color.heatmap = "Reds", direction = 1)
#> Applying a cutoff of  0.05 to the values...

# Take an input of a ligand-receptor pair and show expression in binary
spatialFeaturePlot(cellchat, pairLR.use = "PTHLH_PTH1R", point.size = 1.9, do.binary = TRUE, cutoff = 1.0, enriched.only = F, color.heatmap = "Reds", direction = 1)



