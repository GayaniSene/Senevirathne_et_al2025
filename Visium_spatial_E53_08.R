#Integration of spatial data with single cell data

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

library(spacexr)
library(Seurat)
library(ggplot2)

library(SeuratData)
library(patchwork)
library(dplyr)
library(Matrix)
library(limma)
library(ggstatsplot)
library(ape)
library(CellChat)

install.packages("BiocNeighbors")
BiocManager::install("CellChat")
BiocManager::install("patchwork")

BiocNeighbors

install.packages('BiocManager')
BiocManager::install('BiocNeighbors')

#define the path to the folder that has all the output files from 72_1 batch#1
data_dir <- 'outs_53_08' #dorsal section
data_dir <- 'outs_53_02' #ventral section

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
SpatialFeaturePlot(spatial, features = c("COL9A2", "BMP4"), pt.size = 4)
p1 <- SpatialFeaturePlot(spatial, features = "COL9A2", pt.size.factor = 4)
p2 <- SpatialFeaturePlot(spatial, features = "FGFR3", pt.size.factor = 4, alpha = c(0.1, 1))
p1 + p2


#We can then proceed to run dimensionality reduction and clustering on the RNA expression data, using the same workflow as we use for scRNA-seq analysis.
spatial <- RunPCA(spatial, assay = "SCT", verbose = FALSE)
spatial <- FindNeighbors(spatial, reduction = "pca", dims = 1:30)
spatial <- FindClusters(spatial, verbose = FALSE)
spatial <- RunUMAP(spatial, reduction = "pca", dims = 1:30)


Idents(spatial)


spatial <- RenameIdents(spatial, '0' = 'Mes','1' = 'Mes2','2' = 'Perichondrium')
spatial <- RenameIdents(spatial, '3' = 'Fibo1', '4' ='Chondrocytes', '5' = 'Fibro2')
spatial <- RenameIdents(spatial, '6' = 'MyoProg', '7' = 'Fibo1', '8' ='Chondrocytes')

spatial$celltype <- Idents(spatial)


#save the object as an input for CellChat data.
save(spatial, file = "visium_53_08.RData")


#Umap with only spatial data
p1 <- DimPlot(spatial, reduction = "umap", label = TRUE, pt.size = 3, cols = c('0'='#814827', '1' ='#A47551', '2'= '#60A3D9','3'='#2C5E1A', 
                                                                              '4' = '#055C9D','5'= '#39918C',  '6'='#D43790'))




p2 <- SpatialDimPlot(spatial, label = TRUE, pt.size = 5,label.size = 3, cols = c('0'='#814827', '1' ='#A47551', '2'= '#60A3D9','3'='#2C5E1A', 
                                                                                 '4' = '#055C9D','5'= '#39918C',  '6'='#D43790'))
p1 + p2



plot <- SpatialFeaturePlot(spatial, features = c("Ttr")) + theme(legend.text = element_text(size = 0),
                                                               legend.title = element_text(size = 20), legend.key.size = unit(1, "cm"))
plot <- p1 + p2
jpeg(filename = "UMAP_Spatial_E53.jpg", height = 700, width = 1200, quality = 50)
print(plot)
dev.off()


#highlight genes of interest
SpatialDimPlot(spatial, cells.highlight = CellsByIdentities(object = spatial, idents = c(2, 1, 4, 3,
                                                                                     0, 6)), facet.highlight = TRUE, ncol = 3,pt.size = 3)
de_markers0 <- FindMarkers(spatial, ident.1 = 0)
head(de_markers0, n = 50)

de_markers1 <- FindMarkers(spatial, ident.1 = 1)
head(de_markers1, n = 50)

de_markers2 <- FindMarkers(spatial, ident.1 = 2)
head(de_markers2, n = 50)

de_markers3 <- FindMarkers(spatial, ident.1 = 3)
head(de_markers3, n = 50)

de_markers4 <- FindMarkers(spatial, ident.1 = 4)
head(de_markers4, n = 50)

de_markers5 <- FindMarkers(spatial, ident.1 = 5)
head(de_markers5, n = 50)

de_markers6 <- FindMarkers(spatial, ident.1 = 6)
head(de_markers6, n = 50)

de_markers7 <- FindMarkers(spatial, ident.1 = 7)
head(de_markers7, n = 50)

#defining the cortex with all the cell populations
cortex <- subset(spatial, idents = c(0, 1, 2, 3, 4, 5, 6,7,8))


#Integration with single cell data
allen_reference <- readRDS("pbmc3k_final.rds")

# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k
# cells this speeds up SCTransform dramatically with no loss in performance
library(dplyr)
allen_reference <- SCTransform(allen_reference, ncells = 5000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

# After subsetting, we renormalize cortex
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
# the annotation is stored in the 'subclass' column of object metadata
DimPlot(allen_reference, group.by = "ident", label = TRUE)

anchors <- FindTransferAnchors(reference = allen_reference, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$celltype, prediction.assay = TRUE,
                                  weight.reduction = cortex[["pca"]], dims = 1:30)
cortex[["predictions"]] <- predictions.assay
DefaultAssay(cortex) <- "predictions"

#recognizing the cell types after the integration with single cell data
plot <- SpatialFeaturePlot(cortex, features = c("RestingChondro"), pt.size.factor = 2.2, crop = TRUE)





plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.8), breaks = c(0.0, 0.4, 0.8), low = "lightgray", high = "#055C9D")


plot<-SpatialFeaturePlot(cortex, features = c("Perichondro+Osteoblasts"), pt.size.factor = 5, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0, 0.5), breaks = c(0.0, 0.5), low = "lightgray", high = "#60A3D9")


plot<-SpatialFeaturePlot(cortex, features = c("ChondroProg"), pt.size.factor = 5, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.5), breaks = c(0.0, 0.5), "low" = "lightgray", high = "#003060")


plot<-SpatialFeaturePlot(cortex, features = c("Fibro1"), pt.size.factor = 5, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.10), breaks = c(0.0, 0.10), low = "lightgray", high = "#2C5E1A")


plot<-SpatialFeaturePlot(cortex, features = c("Mes1"), pt.size.factor = 5, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.80), breaks = c(0.0, 0.80), low = "lightgray", high = "#814827")


plot<-SpatialFeaturePlot(cortex, features = c("TransMes"), pt.size.factor = 5,crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.20), breaks = c(0.0, 0.10, 0.20), low = "lightgray", high = "#D2A87E")

plot<-SpatialFeaturePlot(cortex, features = c("Mes2"), pt.size.factor = 5, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.20), breaks = c(0.0, 0.10, 0.20), low = "lightgray", high = "#A47551")

plot<-SpatialFeaturePlot(cortex, features = c("Teno"), pt.size.factor = 5, ncol = 2, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.20), breaks = c(0.0, 0.10, 0.20), low = "lightgray", high = "#F2C5E0")




plot<-SpatialFeaturePlot(cortex, features = c("Fibro2"), pt.size.factor = 5, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.20), breaks = c(0.0, 0.10, 0.20), low = "lightgray", high = "#39918C")


plot<-SpatialFeaturePlot(cortex, features = c("MyoProg"), pt.size.factor = 5, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.20), breaks = c(0.0, 0.10, 0.20), low = "lightgray", high = "#D43790")

plot<-SpatialFeaturePlot(cortex, features = c("Fibro3"), pt.size.factor = 5, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.10), breaks = c(0.0, 0.10), low = "lightgray", high = "#B2D2A4")

plot<-SpatialFeaturePlot(cortex, features = c("RestingChondro2"), pt.size.factor = 5, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.30), breaks = c(0.0, 0.10, 0.30), low = "lightgray", high = "#055C9D")

plot<-SpatialFeaturePlot(cortex, features = c("Perimysium"), pt.size.factor = 5, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.20), breaks = c(0.0, 0.20), low = "lightgray", high = "#FF0080")


plot<-SpatialFeaturePlot(cortex, features = c("ProxMes"), pt.size.factor = 5, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.20), breaks = c(0.0, 0.10, 0.20), low = "lightgray", high = "#D0B49F")

plot<-SpatialFeaturePlot(cortex, features = c("MyoProg+Pax3"), pt.size.factor = 5, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.20), breaks = c(0.0, 0.10, 0.20), low = "lightgray", high = "#EC8FD0")


plot<-SpatialFeaturePlot(cortex, features = c("SmoothMuscleProg"), pt.size.factor = 5, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.20), breaks = c(0.0, 0.10, 0.20), low = "lightgray", high = "#FB8DA0")

plot<-SpatialFeaturePlot(cortex, features = c("EndothelialCells"), pt.size.factor = 5, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.20), breaks = c(0.0, 0.10, 0.20), low = "lightgray", high = "#973B10")

plot<-SpatialFeaturePlot(cortex, features = c("Schwann"), pt.size.factor = 5, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.20), breaks = c(0.0, 0.10, 0.20), low = "lightgray", high = "#FFDD64")

plot<-SpatialFeaturePlot(cortex, features = c("ImmuneCells"), pt.size.factor = 5, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.20), breaks = c(0.0, 0.10, 0.20), low = "lightgray", high = "#603F8B")








