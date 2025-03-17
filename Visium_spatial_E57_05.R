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

install.packages("ggstatsplot")

install.packages('BiocManager')
BiocManager::install('limma')

#define the path to the folder that has all the output files from 72_1 batch#1
data_dir <- 'out_57_05'
data_dir <- 'outs_57_09'
data_dir <- 'outs_57_01_new'
data_dir <- 'outs_57_02_new'
data_dir <- 'outs_57_03_new'
data_dir <- 'outs_57_04_new'

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
p2 <- SpatialFeaturePlot(spatial, features = "RUNX2", pt.size.factor = 3, alpha = c(0.1, 1))
p2 <- SpatialFeaturePlot(spatial, features = "UCMA", pt.size.factor = 3, alpha = c(0.1, 1))
p2 <- SpatialFeaturePlot(spatial, features = "SP7", pt.size.factor = 3, alpha = c(0.1, 1))

p1 + p2


#We can then proceed to run dimensionality reduction and clustering on the RNA expression data, using the same workflow as we use for scRNA-seq analysis.
spatial <- RunPCA(spatial, assay = "SCT", verbose = FALSE)
spatial <- FindNeighbors(spatial, reduction = "pca", dims = 1:30)
spatial <- FindClusters(spatial, verbose = FALSE)
spatial <- RunUMAP(spatial, reduction = "pca", dims = 1:30)


pbmc <- RenameIdents(spatial, '0' = 'PrimaryErythro','1' = 'Fibro1','2' = 'MyoProg')
pbmc <- RenameIdents(spatial, '3' = 'Chondrocytes', '4' ='Fibro2', '5' = 'Muscles')
pbmc <- RenameIdents(spatial, '6' = 'Mes', '7'='periderm', '8' = 'osteoblasts', '9'='MyoC', '10'='schwann', '11'='perimysium', '12'='Mes2')

pbmc$celltype <- Idents(pbmc)


pbmc <- RenameIdents(spatial, '0' = 'PrimaryErythro','1' = 'Fibro1','2' = 'MyoProg')
pbmc <- RenameIdents(spatial, '3' = 'Chondrocytes', '4' ='Fibro2', '5' = 'Muscles')
pbmc <- RenameIdents(spatial, '6' = 'Mes', '7'='periderm', '8' = 'osteoblasts')


#Umap with only spatial data
p1 <- DimPlot(spatial, reduction = "umap", label = TRUE, pt.size = 3, cols = c('0'='#FD7F20', '1' ='#2C5E1A', '2'= '#D43790','3'='#055C9D', 
                                                                               '4' = '#39918C','5'= '#EC8FD0',  '6'='#814827', '7'='#D2A87E', 
                                                                                '8'='#60A3D9', '9'= '#FB4570', '10'='#FFDD64', '11'='#FF0080', '12'='#A47551'))




p2 <- SpatialDimPlot(spatial, label = TRUE, pt.size = 3,label.size = 3, cols = c('0'='#FD7F20', '1' ='#2C5E1A', '2'= '#D43790','3'='#055C9D', 
                                                                                 '4' = '#39918C','5'= '#EC8FD0',  '6'='#814827', '7'='#D2A87E', 
                                                                                 '8'='#60A3D9', '9'= '#FB4570', '10'='#FFDD64', '11'='#FF0080', '12'='#A47551'))
p1 + p2



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
head(de_markers7, n = 20)

de_markers8 <- FindMarkers(spatial, ident.1 = 8)
head(de_markers8, n = 20)

de_markers9 <- FindMarkers(spatial, ident.1 = 9)
head(de_markers9, n = 20)

de_markers10 <- FindMarkers(spatial, ident.1 = 10)
head(de_markers10, n = 20)

de_markers11 <- FindMarkers(spatial, ident.1 = 11)
head(de_markers11, n = 20)

de_markers12 <- FindMarkers(spatial, ident.1 = 12)
head(de_markers12, n = 20)

#defining the cortex with all the cell populations
cortex <- subset(spatial, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11))
cortex <- subset(spatial, idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8))

#Integration with single cell data
allen_reference <- readRDS("pbmc3k_final.rds")

# note that setting ncells=3000 normalizes the full dataset but learns noise models on 3k
# cells this speeds up SCTransform dramatically with no loss in performance
library(dplyr)
allen_reference <- SCTransform(allen_reference, ncells = 6000, verbose = FALSE) %>%
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


plot1 <- SpatialFeaturePlot(cortex, features = c("RestingChondro"), pt.size.factor = 2.8, crop = TRUE)
plot2 <- SpatialFeaturePlot(cortex, features = c("Myoprog+Pax7"), pt.size.factor = 2.8, crop = TRUE)
plot3 <- SpatialFeaturePlot(cortex, features = c("Perichondrium"), pt.size.factor = 2.8, crop = TRUE)
plot4 <- SpatialFeaturePlot(cortex, features = c("Proliferating"), pt.size.factor = 2.8, crop = TRUE)
plot5 <- SpatialFeaturePlot(cortex, features = c("Mes2"), pt.size.factor = 2.8, crop = TRUE)
plot6 <- SpatialFeaturePlot(cortex, features = c("Mes1"), pt.size.factor = 2.8, crop = TRUE)
plot8 <- SpatialFeaturePlot(cortex, features = c("Fibro1"), pt.size.factor = 2.8, crop = TRUE)
plot9 <- SpatialFeaturePlot(cortex, features = c("Perimysium"), pt.size.factor = 2.8, crop = TRUE)
plot10 <- SpatialFeaturePlot(cortex, features = c("PrimaryErythro"), pt.size.factor = 2.8, crop = TRUE)
plot11 <- SpatialFeaturePlot(cortex, features = c("Erythro"), pt.size.factor = 2.8, crop = TRUE)
plot12 <- SpatialFeaturePlot(cortex, features = c("Proliferating2"), pt.size.factor = 2.8, crop = TRUE)
plot13 <- SpatialFeaturePlot(cortex, features = c("Proliferating3"), pt.size.factor = 2.8, crop = TRUE)
plot14 <- SpatialFeaturePlot(cortex, features = c("TransMes"), pt.size.factor = 2.8, crop = TRUE)
plot15 <- SpatialFeaturePlot(cortex, features = c("VenousEndo"), pt.size.factor = 2.8, crop = TRUE)
plot16 <- SpatialFeaturePlot(cortex, features = c("ArterialEndo"), pt.size.factor = 2.8, crop = TRUE)
plot17 <- SpatialFeaturePlot(cortex, features = c("MyoC"), pt.size.factor = 2.8, crop = TRUE)
plot18 <- SpatialFeaturePlot(cortex, features = c("Fibro2"), pt.size.factor = 2.8, crop = TRUE)
plot19 <- SpatialFeaturePlot(cortex, features = c("Macrophages"), pt.size.factor = 2.8, crop = TRUE)
plot20 <- SpatialFeaturePlot(cortex, features = c("Schwann"), pt.size.factor = 2.8, crop = TRUE)
plot7 <- SpatialFeaturePlot(cortex, features = c("perichondrium_Osteoblast"), pt.size.factor = 2.8, crop = TRUE)



plot <- SpatialFeaturePlot(cortex, features = c("RestingChondro"), pt.size.factor = 3, crop = TRUE,  alpha = c(0.7,1))
plot + scale_fill_continuous(limits = c(0.0, 0.025), breaks = c(0.0, 0.025), low = "lightgray", high = "#055C9D")

plot <- SpatialFeaturePlot(cortex, features = c("RestingChondro"), pt.size.factor = 3, crop = TRUE,alpha = c(0.9,1))
plot + scale_fill_continuous(limits = c(0.0, 0.025), breaks = c(0.0, 0.025), low = "darkgray", high = "#055C9D")


alpha = c(0.1, 1)
plot<-SpatialFeaturePlot(cortex, features = c("Myoprog+Pax7"), pt.size.factor = 3, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.1), breaks = c(0.0, 0.10), low = "lightgray", high = "#D43790")


plot<-SpatialFeaturePlot(cortex, features = c("Perichondrium"), pt.size.factor = 3, crop = TRUE, alpha = c(0.9,1))
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.8), breaks = c(0.0, 0.40, 0.80), low = "darkgray", high = "#60A3D9")


plot<-SpatialFeaturePlot(cortex, features = c("Proliferating"), pt.size.factor = 3, crop = TRUE, alpha = c(0.9, 1))
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.8), breaks = c(0.0,0.45), low = "darkgray", high = "#0E86D4")


plot<-SpatialFeaturePlot(cortex, features = c("Mes2"), pt.size.factor = 3, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.20), breaks = c(0.0, 0.10, 0.20), low = "lightgray", high = "#A47551")

plot<-SpatialFeaturePlot(cortex, features = c("Mes1"), pt.size.factor = 3,crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.10), breaks = c(0.0, 0.10, 0.20), low = "lightgray", high = "#814827")


#Not found
plot<-SpatialFeaturePlot(cortex, features = c("perichondrium_Osteoblast"), pt.size.factor = 3, crop = TRUE,alpha = c(0.7, 1))
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.90), breaks = c(0.0, 0.10, 0.20), low = "lightgray", high = "#68BBE3")


plot<-SpatialFeaturePlot(cortex, features = c("Fibro1"), pt.size.factor = 3, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.10), breaks = c(0.0, 0.10, 0.20), low = "lightgray", high = "#2C5E1A")

plot<-SpatialFeaturePlot(cortex, features = c("Perimysium"), pt.size.factor = 3, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.10), breaks = c(0.0, 0.10, 0.20), low = "lightgray", high = "#FF0080")

plot<-SpatialFeaturePlot(cortex, features = c("PrimaryErythro"), pt.size.factor = 3, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,.60), breaks = c(0.0, 0.30, 0.60), low = "lightgray", high = "#F5631A")

plot<-SpatialFeaturePlot(cortex, features = c("Erythro"), pt.size.factor = 3, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.90), breaks = c(0.0, 0.90), low = "lightgray", high = "#FD7F20")

plot<-SpatialFeaturePlot(cortex, features = c("Proliferating2"), pt.size.factor = 3, crop = TRUE, alpha = c(0.7, 1))
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.05), breaks = c(0.0, 0.03), low = "lightgray", high = "#0E86D4")


plot<-SpatialFeaturePlot(cortex, features = c("Proliferating3"), pt.size.factor = 3, crop = TRUE, alpha = c(0.9,1))
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.5), breaks = c(0.0, 0.5), low = "darkgray", high = "#0E86D4")


plot<-SpatialFeaturePlot(cortex, features = c("TransMes"), pt.size.factor = 3, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.10), breaks = c(0.0, 0.10, 0.20), low = "lightgray", high = "#D2A87E")

plot<-SpatialFeaturePlot(cortex, features = c("VenousEndo"), pt.size.factor = 3, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.10), breaks = c(0.0, 0.10, 0.20), low = "lightgray", high = "#D85E00")


plot<-SpatialFeaturePlot(cortex, features = c("MyoC"), pt.size.factor = 3, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.90), breaks = c(0.0, 0.45, 0.90), low = "lightgray", high = "#FB4570")

plot<-SpatialFeaturePlot(cortex, features = c("ArterialEndo"), pt.size.factor = 3, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.20), breaks = c(0.0, 0.10, 0.20), low = "lightgray", high = "#E55B13")

plot<-SpatialFeaturePlot(cortex, features = c("Fibro2"), pt.size.factor = 3, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.10), breaks = c(0.0, 0.10, 0.20), low = "lightgray", high = "#39918C")

plot<-SpatialFeaturePlot(cortex, features = c("Macrophages"), pt.size.factor = 3, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.90), breaks = c(0.0, 0.90), low = "lightgray", high = "#603F8B")

plot<-SpatialFeaturePlot(cortex, features = c("Schwann"), pt.size.factor = 3, crop = TRUE)
plot + ggplot2::scale_fill_continuous(limits = c(0.0,0.90), breaks = c(0.0, 0.90), low = "lightgray", high = "#FFDD64")







