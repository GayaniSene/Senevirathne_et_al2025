#R code for the E53_IL_new_Multiomics



#install the packages
install.packages('Seurat')
library(Seurat)
setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
install.packages(c("BPCells", "presto", "glmGamPoi"))
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
install.packages('Signac')
install.packages("devtools")
devtools::install_github("chromVAR", build_vignettes = FALSE)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")
remotes::install_version("Matrix", version = "1.6-1")
packageVersion("gfortran")
remotes::install_version("Matrix", version = "1.6-1.1")
remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
install.packages("Seurat")
install.packages("ggplot2")
install.packages("patchwork")
install.packages("dplyr")
install.packages("limma")
install.packages("ggstatsplot")
install.packages("slam")
install.packages("Signac")
install.packages("GenomeInfoDb")
install.packages("EnsDb.Hsapiens.v75")
install.packages("EnsDb.Hsapiens.v86")
install.packages("hdf5r")
set.seed(1234)
install.packages("biovizBase")
install.packages("stringr")
require(devtools)
install.packages("tidyverse")
install.packages("BSgenome.Hsapiens.UCSC.hg19")
install.packages("chromVAR")
install.packages("JASPAR2020")
install.packages("TFBSTools")
install.packages("motifmatchr")
install.packages("BSgenome.Hsapiens.1000genomes.hs37d5")
install.packages("presto")
install.packages("ggseqlogo")
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(limma)
library(ggstatsplot)
library(slam)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(hdf5r)
set.seed(1234)
library(biovizBase)
library(stringr)
require(devtools)
library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg19)
library(chromVAR)
library(TFBSTools)
library(presto)
library("ggseqlogo")
library(CellChat)

library('ggalluvial')
library(CellChat)
library(patchwork)


# the 10x hdf5 file contains both data types. 
inputdata.10x <- Read10X_h5("filtered_feature_bc_matrix.h5")

# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

# Create Seurat object
IL <- CreateSeuratObject(counts = rna_counts)
IL[["percent.mt"]] <- PercentageFeatureSet(IL, pattern = "^MT-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"


chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = "atac_fragments.tsv.gz",
  min.cells = 10,
  annotation = annotations
)

IL[["ATAC"]] <- chrom_assay

#We perform basic QC based on the number of detected molecules for each modality as well as mitochondrial percentage.

VlnPlot(IL, features = c("nCount_ATAC", "nCount_RNA", "percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()

#If you want to subset the dataset based on the QC values,
IL <- subset(
  x = IL,
  subset = nCount_ATAC < XX &
    nCount_ATAC > XX &
    nCount_RNA < XX &
    nCount_RNA > XX &
    percent.mt < XX
)

# RNA analysis
DefaultAssay(IL) <- "RNA"
IL <- SCTransform(IL, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(IL) <- "ATAC"
IL <- RunTFIDF(IL)
IL <- FindTopFeatures(IL, min.cutoff = 'q0')
IL <- RunSVD(IL)
IL <- RunUMAP(IL, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

#We calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. We use this graph for UMAP visualization and clustering

IL <- FindMultiModalNeighbors(IL, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
IL <- RunUMAP(IL, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
IL <- FindClusters(IL, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

#clustering based on gene expression, ATAC-seq, or WNN analysis.
p1 <- DimPlot(IL, reduction = "umap.rna", label = FALSE, label.size = 3.5, pt.size = 1.4, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(IL, reduction = "umap.atac",  label = FALSE, label.size = 3.5, pt.size = 1.4, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(IL, reduction = "wnn.umap",  label = TRUE, label.size = 3.5, pt.size = 2, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

p3

#Renaming the clusters

IL <- RenameIdents(IL, '0' = 'ChondroProg','1' = 'Fibro1','2' = 'Mes1')
IL <- RenameIdents(IL, '3' = 'TransMes', '4' ='Mes2', '5' = 'Teno')
IL <- RenameIdents(IL, '6' = 'Fibro2','7' = 'MyoProg', '8' = 'Perichondro+Osteoblasts')
IL <- RenameIdents(IL, '9' = 'RestingChondro', '10' = 'Fibro3', '11' = 'RestingChondro2')
IL <- RenameIdents(IL, '12' = 'Perimysium', '13' = 'ProxMes', '14' = 'MyoProg+Pax3')
IL <- RenameIdents(IL, '15' = 'SmoothMuscleProg', '16' = 'EndothelialCells', '17' = 'Schwann')
IL <- RenameIdents(IL, '18' = 'ImmuneCells')
IL$celltype <- Idents(IL)

save(IL, file = "E53.RData")


#clustering based on gene expression, ATAC-seq, or WNN analysis and colored according to cell type.
p1 <- DimPlot(IL, reduction = "umap.rna", label = FALSE, label.size = 3.5, pt.size = 2.5, repel = TRUE,cols = c('RestingChondro2'='#055C9D', 'Fibro1'= '#2C5E1A', 'Mes1'='#814827', 'TransMes'= '#D2A87E', 'Mes2' = '#A47551', 'Teno'='#F2C5E0',
                                                                                                                  'Fibro2'= '#39918C', 'MyoProg'= '#D43790' , 'Perichondro+Osteoblasts'= '#68BBE3', 'ChondroProg'= '#003060', 'Fibro3' = '#B2D2A4', 'RestingChondro'='#055C9D',
                                                                                                                  'Perimysium' = '#FF0080', 'ProxMes'= '#D0B49F', 'MyoProg+Pax3'= '#EC8FD0', 'SmoothMuscleProg'= '#FB8DA0', 'EndothelialCells'= '#973B10',
                                                                                                                  'Schwann'= '#FFDD64','ImmuneCells'='#603F8B')) + ggtitle("RNA")
p2 <- DimPlot(IL, reduction = "umap.atac",  label = FALSE, label.size = 3.5, pt.size = 2.5, repel = TRUE, cols = c('RestingChondro2'='#055C9D', 'Fibro1'= '#2C5E1A', 'Mes1'='#814827', 'TransMes'= '#D2A87E', 'Mes2' = '#A47551', 'Teno'='#F2C5E0',
                                                                                                                     'Fibro2'= '#39918C', 'MyoProg'= '#D43790' , 'Perichondro+Osteoblasts'= '#68BBE3', 'ChondroProg'= '#003060', 'Fibro3' = '#B2D2A4', 'RestingChondro'='#055C9D',
                                                                                                                     'Perimysium' = '#FF0080', 'ProxMes'= '#D0B49F', 'MyoProg+Pax3'= '#EC8FD0', 'SmoothMuscleProg'= '#FB8DA0', 'EndothelialCells'= '#973B10',
                                                                                                                     'Schwann'= '#FFDD64','ImmuneCells'='#603F8B')) + ggtitle("ATAC")
p3 <- DimPlot(IL, reduction = "wnn.umap",  label = TRUE, label.size = 3.5, pt.size = 2, repel = TRUE) + ggtitle("WNN")

p3 <- DimPlot(IL, reduction = "wnn.umap",  label = TRUE, label.size = 4, pt.size = 3.5, repel = TRUE, cols = c('RestingChondro2'='#055C9D', 'Fibro1'= '#2C5E1A', 'Mes1'='#814827', 'TransMes'= '#D2A87E', 'Mes2' = '#A47551', 'Teno'='#F2C5E0',
                                                                                                                 'Fibro2'= '#39918C', 'MyoProg'= '#D43790' , 'Perichondro+Osteoblasts'= '#68BBE3', 'ChondroProg'= '#003060', 'Fibro3' = '#B2D2A4', 'RestingChondro'='#055C9D',
                                                                                                                 'Perimysium' = '#FF0080', 'ProxMes'= '#D0B49F', 'MyoProg+Pax3'= '#EC8FD0', 'SmoothMuscleProg'= '#FB8DA0', 'EndothelialCells'= '#973B10',
                                                                                                                 'Schwann'= '#FFDD64','ImmuneCells'='#603F8B')) + ggtitle("WNN")

p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
p3

DefaultAssay(IL) <- "SCT"
# we can also identify alternative protein and RNA markers for this cluster through
# differential expression

rna_markers0 <- FindMarkers(IL, ident.1 = 0, assay = "SCT")
head(rna_markers0)
head(rna_markers0, n = 50)

rna_markers1 <- FindMarkers(IL, ident.1 = 1, assay = "SCT")
head(rna_markers1)
head(rna_markers1, n = 50)

rna_markers2 <- FindMarkers(IL, ident.1 = 2, assay = "SCT")
head(rna_markers2)
head(rna_markers2, n = 50)

rna_markers3 <- FindMarkers(IL, ident.1 = 3, assay = "SCT")
head(rna_markers3)
head(rna_markers3, n = 50)

rna_markers4 <- FindMarkers(IL, ident.1 = 4, assay = "SCT")
head(rna_markers4)
head(rna_markers4, n = 50)

rna_markers5 <- FindMarkers(IL, ident.1 = 5, assay = "SCT")
head(rna_markers5)
head(rna_markers5, n = 50)

rna_markers6 <- FindMarkers(IL, ident.1 = 6, assay = "SCT")
head(rna_markers6)
head(rna_markers6, n = 50)

rna_markers7 <- FindMarkers(IL, ident.1 = 7, assay = "SCT")
head(rna_markers7)
head(rna_markers7, n = 50)

rna_markers8 <- FindMarkers(IL, ident.1 = 8, assay = "SCT")
head(rna_markers8)
head(rna_markers8, n = 50)

rna_markers9 <- FindMarkers(IL, ident.1 = 9, assay = "SCT")
head(rna_markers9)
head(rna_markers9, n = 50)

rna_markers10 <- FindMarkers(IL, ident.1 = 10, assay = "SCT")
head(rna_markers10)
head(rna_markers10, n = 50)

rna_markers11 <- FindMarkers(IL, ident.1 = 11, assay = "SCT")
head(rna_markers11)
head(rna_markers11, n = 50)

rna_markers12 <- FindMarkers(IL, ident.1 = 12, assay = "SCT")
head(rna_markers12)
head(rna_markers12, n = 50)

rna_markers13 <- FindMarkers(IL, ident.1 = 13, assay = "SCT")
head(rna_markers13)
head(rna_markers13, n = 50)

rna_markers14 <- FindMarkers(IL, ident.1 = 14, assay = "SCT")
head(rna_markers14)
head(rna_markers14, n = 50)

rna_markers15 <- FindMarkers(IL, ident.1 = 15, assay = "SCT")
head(rna_markers15)
head(rna_markers15, n = 50)

rna_markers16 <- FindMarkers(IL, ident.1 = 16, assay = "SCT")
head(rna_markers16)
head(rna_markers16, n = 50)

rna_markers17 <- FindMarkers(IL, ident.1 = 17, assay = "SCT")
head(rna_markers17)
head(rna_markers17, n = 50)

rna_markers18 <- FindMarkers(IL, ident.1 = 18, assay = "SCT")
head(rna_markers18)
head(rna_markers18, n = 50)

rna_markers19 <- FindMarkers(IL, ident.1 = 19, assay = "SCT")
head(rna_markers19)
head(rna_markers19, n = 20)

#save the RDS file
saveRDS(IL, file = "IL3k_final.rds")

#marker gene plots in Extended Data Fig3. 
DefaultAssay(IL) <- "SCT"
features <- c("MATN4","CNMD","UCMA","COL9A2","MATN3","COMP","SNORC","SHOX","DLX5", "RUNX2", "FGFBP2","IHH","THBS2","COL2A1","COL1A2","COL10A1", "PAX3", "PAX7", "TNMD", "IGFBP7", "MPZ", "SP7", "SOX9", "VEGFA", "VEGFB", "ZNF521", "PTH1R","FOXD3")

Features <-("UCMA")
DotPlot(IL, features = features) + RotatedAxis()
FeaturePlot(IL, features = "COL9A2") + RotatedAxis()

FeaturePlot(IL, features = "sct_WWP2", reduction = 'wnn.umap', pt.size = 1.5)
FeaturePlot(IL, features = "sct_SOX9", reduction = 'wnn.umap', pt.size = 1.5)
FeaturePlot(IL, features = "sct_UCMA", reduction = 'wnn.umap', pt.size = 1.5)
FeaturePlot(IL, features = "sct_COL2A1", reduction = 'wnn.umap', pt.size = 1.5)
FeaturePlot(IL, features = "sct_COL9A2", reduction = 'wnn.umap', pt.size = 1.5)
FeaturePlot(IL, features = "sct_PTH1R", reduction = 'wnn.umap', pt.size = 1.5)
FeaturePlot(IL, features = "sct_IHH", reduction = 'wnn.umap', pt.size = 1.5)
FeaturePlot(IL, features = "sct_COL10A1", reduction = 'wnn.umap', pt.size = 1.5)


FeaturePlot(IL, features = "sct_THBS2", reduction = 'wnn.umap', pt.size = 1.5)

FeaturePlot(IL, features = "sct_RUNX2", reduction = 'wnn.umap', pt.size = 1.5)
FeaturePlot(IL, features = "sct_SP7", reduction = 'wnn.umap', pt.size = 1.5)
FeaturePlot(IL, features = "sct_DLX5", reduction = 'wnn.umap', pt.size = 1.5)

FeaturePlot(IL, features = "sct_VEGFA", reduction = 'wnn.umap', pt.size = 1.5)
FeaturePlot(IL, features = "sct_VEGFB", reduction = 'wnn.umap', pt.size = 1.5)
FeaturePlot(IL, features = "sct_VEGFC", reduction = 'wnn.umap', pt.size = 1.5)
FeaturePlot(IL, features = "sct_FOXP1", reduction = 'wnn.umap', pt.size = 1.5)
FeaturePlot(IL, features = "sct_FOXP2", reduction = 'wnn.umap', pt.size = 1.5)


