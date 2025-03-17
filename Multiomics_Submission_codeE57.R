#R code for the E57_IL_Multiomics

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
IL_57 <- CreateSeuratObject(counts = rna_counts)
IL_57[["percent.mt"]] <- PercentageFeatureSet(IL_57, pattern = "^MT-")

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

IL_57[["ATAC"]] <- chrom_assay


#We perform basic QC based on the number of detected molecules for each modality as well as mitochondrial percentage.

VlnPlot(IL_57, features = c("nCount_ATAC", "nCount_RNA", "percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()

#If you want to subset the dataset based on the QC values,
IL_57 <- subset(
  x = IL_57,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)

# RNA analysis
DefaultAssay(IL_57) <- "RNA"
IL_57 <- SCTransform(IL_57, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(IL_57) <- "ATAC"
IL_57 <- RunTFIDF(IL_57)
IL_57 <- FindTopFeatures(IL_57, min.cutoff = 'q0')
IL_57 <- RunSVD(IL_57)
IL_57 <- RunUMAP(IL_57, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

#We calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. We use this graph for UMAP visualization and clustering

IL_57 <- FindMultiModalNeighbors(IL_57, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
IL_57 <- RunUMAP(IL_57, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
IL_57 <- FindClusters(IL_57, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

#clustering based on gene expression, ATAC-seq, or WNN analysis.
p1 <- DimPlot(IL_57, reduction = "umap.rna", label = FALSE, label.size = 3.5, pt.size = 1.4, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(IL_57, reduction = "umap.atac",  label = FALSE, label.size = 3.5, pt.size = 1.4, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(IL_57, reduction = "wnn.umap",  label = TRUE, label.size = 3.5, pt.size = 2, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

p3
jpeg(filename = "UMAP_E57_IL.jpg", quality = 500)
print(p3)
dev.off()

p3
tiff('test.tiff', units="in", width=10, height=8, res=300, compression = 'lzw')
print(p3)
dev.off()


#cell types are names as idents before labelling them (this worked for me)
IL_57$celltype = Idents(IL_57)

# we can also identify alternative protein and RNA markers for this cluster through
# differential expression
rna_markers0 <- FindMarkers(IL_57, ident.1 = 0, assay = "SCT")
head(rna_markers0)
head(rna_markers0, n = 50)


rna_markers1 <- FindMarkers(IL_57, ident.1 = 1, assay = "SCT")
head(rna_markers1)
head(rna_markers1, n = 50)

rna_markers2 <- FindMarkers(IL_57, ident.1 = 2, assay = "SCT")
head(rna_markers2)
head(rna_markers2, n = 50)

rna_markers3 <- FindMarkers(IL_57, ident.1 = 3, assay = "SCT")
head(rna_markers3)
head(rna_markers3, n = 50)

rna_markers4 <- FindMarkers(IL_57, ident.1 = 4, assay = "SCT")
head(rna_markers4)
head(rna_markers4, n = 50)

rna_markers5 <- FindMarkers(IL_57, ident.1 = 5, assay = "SCT")
head(rna_markers5)
head(rna_markers5, n = 50)

rna_markers6 <- FindMarkers(IL_57, ident.1 = 6, assay = "SCT")
head(rna_markers6)
head(rna_markers6, n = 50)

rna_markers7 <- FindMarkers(IL_57, ident.1 = 7, assay = "SCT")
head(rna_markers7)
head(rna_markers7, n = 50)

rna_markers8 <- FindMarkers(IL_57, ident.1 = 8, assay = "SCT")
head(rna_markers8)
head(rna_markers8, n = 50)

rna_markers9 <- FindMarkers(IL_57, ident.1 = 9, assay = "SCT")
head(rna_markers9)
head(rna_markers9, n = 50)

rna_markers10 <- FindMarkers(IL_57, ident.1 = 10, assay = "SCT")
head(rna_markers10)
head(rna_markers10, n = 50)

rna_markers11 <- FindMarkers(IL_57, ident.1 = 11, assay = "SCT")
head(rna_markers11)
head(rna_markers11, n = 50)

rna_markers12 <- FindMarkers(IL_57, ident.1 = 12, assay = "SCT")
head(rna_markers12)
head(rna_markers12, n = 50)

rna_markers13 <- FindMarkers(IL_57, ident.1 = 13, assay = "SCT")
head(rna_markers13)
head(rna_markers13, n = 50)

rna_markers14 <- FindMarkers(IL_57, ident.1 = 14, assay = "SCT")
head(rna_markers14)
head(rna_markers14, n = 50)

rna_markers15 <- FindMarkers(IL_57, ident.1 = 15, assay = "SCT")
head(rna_markers15)
head(rna_markers15, n = 50)

rna_markers16 <- FindMarkers(IL_57, ident.1 = 16, assay = "SCT")
head(rna_markers16)
head(rna_markers16, n = 50)

rna_markers17 <- FindMarkers(IL_57, ident.1 = 17, assay = "SCT")
head(rna_markers17)
head(rna_markers17, n = 50)

rna_markers18 <- FindMarkers(IL_57, ident.1 = 18, assay = "SCT")
head(rna_markers18)
head(rna_markers18, n = 50)

rna_markers19 <- FindMarkers(IL_57, ident.1 = 19, assay = "SCT")
head(rna_markers19)
head(rna_markers19, n = 20)



#Renaming the clusters
IL_57 <- RenameIdents(IL_57, '0' = 'RestingChondro')
IL_57 <- RenameIdents(IL_57, '1' = 'Myoprog+Pax7') 
IL_57 <- RenameIdents(IL_57, '2' = 'Perichondrium') 
IL_57 <- RenameIdents(IL_57, '3' = 'Proliferating')
IL_57 <- RenameIdents(IL_57, '4' ='Mes2')
IL_57 <- RenameIdents(IL_57, '5' = 'Mes1') 
IL_57 <- RenameIdents(IL_57, '6' = 'Osteoblast') 
IL_57 <- RenameIdents(IL_57, '7' = 'Fibro1')
IL_57 <- RenameIdents(IL_57, '8' = 'Perimysium') 
IL_57 <- RenameIdents(IL_57, '9' = 'PrimaryErythro') 
IL_57 <- RenameIdents(IL_57, '10' = 'Erythro') 
IL_57 <- RenameIdents(IL_57, '11' = 'Proliferating2')
IL_57 <- RenameIdents(IL_57, '12' = 'Proliferating3') 
IL_57 <- RenameIdents(IL_57, '13' = 'TransMes')
IL_57 <- RenameIdents(IL_57, '14' = 'VenousEndo') 
IL_57 <- RenameIdents(IL_57, '15' = 'MyoC', '16' = 'ArterialEndo', '17'= 'Fibro2', '18' = 'Macrophages', '19'='Schwann')


IL_57$celltype <- Idents(IL_57)


p1 <- DimPlot(IL_57, reduction = "umap.rna", label = FALSE, label.size = 3.5, pt.size = 1.4, repel = TRUE, cols = c('RestingChondro' = '#055C9D', 'Myoprog+Pax7' = '#D43790', 'Perichondrium' = '#60A3D9', 'Proliferating' = '#0E86D4', 'Mes2' = '#A47551',
                                                                                                                   'Mes1' = '#814827', 'Osteoblast' = '#BA0F30', 'Fibro1' = '#2C5E1A','Perimysium' = '#FF0080','PrimaryErythro' = '#F5631A',
                                                                                                                   'Erythro' =  '#FD7F20','Proliferating2' = '#0E86D4', 'Proliferating3' =  '#0E86D4','TransMes' = '#D2A87E', 'VenousEndo' =  '#D85E00',
                                                                                                                   'MyoC'= '#FB4570','ArterialEndo' = '#E55B13', 'Fibro2' = '#39918C','Macrophages' = '#603F8B','Schwann' = '#FFDD64')) + ggtitle("RNA")
p2 <- DimPlot(IL_57, reduction = "umap.atac",  label = FALSE, label.size = 3.5, pt.size = 1.4, repel = TRUE, cols = c('RestingChondro' = '#055C9D', 'Myoprog+Pax7' = '#D43790', 'Perichondrium' = '#60A3D9', 'Proliferating' = '#0E86D4', 'Mes2' = '#A47551',
                                                                                                                     'Mes1' = '#814827', 'Osteoblast' = '#BA0F30', 'Fibro1' = '#2C5E1A','Perimysium' = '#FF0080','PrimaryErythro' = '#F5631A',
                                                                                                                     'Erythro' =  '#FD7F20','Proliferating2' = '#0E86D4', 'Proliferating3' =  '#0E86D4','TransMes' = '#D2A87E', 'VenousEndo' =  '#D85E00',
                                                                                                                     'MyoC'= '#FB4570','ArterialEndo' = '#E55B13', 'Fibro2' = '#39918C','Macrophages' = '#603F8B','Schwann' = '#FFDD64')) + ggtitle("ATAC")
p3 <- DimPlot(IL_57, reduction = "wnn.umap",  label = TRUE, label.size = 3.5, pt.size = 2, repel = TRUE, cols = c('RestingChondro' = '#055C9D', 'Myoprog+Pax7' = '#D43790', 'Perichondrium' = '#60A3D9', 'Proliferating' = '#0E86D4', 'Mes2' = '#A47551',
                                                                                                                 'Mes1' = '#814827', 'Osteoblast' = '#BA0F30', 'Fibro1' = '#2C5E1A','Perimysium' = '#FF0080','PrimaryErythro' = '#F5631A',
                                                                                                                 'Erythro' =  '#FD7F20','Proliferating2' = '#0E86D4', 'Proliferating3' =  '#0E86D4','TransMes' = '#D2A87E', 'VenousEndo' =  '#D85E00',
                                                                                                                 'MyoC'= '#FB4570','ArterialEndo' = '#E55B13', 'Fibro2' = '#39918C','Macrophages' = '#603F8B','Schwann' = '#FFDD64')) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))


#save the RDS file
saveRDS(IL_57, file = "IL_573k_final.rds")


#marker gene plots in Extended Data Fig3. 

FeaturePlot(IL_57, features = "sct_NCL", reduction = 'wnn.umap', pt.size = 1)
FeaturePlot(IL_57, features = "sct_SDC2", reduction = 'wnn.umap', pt.size = 1)
FeaturePlot(IL_57, features = "sct_", reduction = 'wnn.umap', pt.size = 1)


DefaultAssay(IL_57) <- "SCT"
features <- c("MATN4","CNMD","UCMA","COL9A2","MATN3","COMP","SNORC","SHOX","DLX5", "RUNX2", "FGFBP2","IHH","THBS2","COL2A1","COL1A2","COL10A1", "PAX3", "PAX7", "TNMD", "IGFBP7", "MPZ", "SP7", "SOX9", "VEGFA", "VEGFB", "ZNF521", "PTH1R","FOXD3")

DotPlot(IL_57, features = features) + RotatedAxis()

              
features <- c("FRZB","BMP6","UCMA","CNMD","DLK1","MYF5","PAX7","PDLIM3","VGLL3","TSPAN12","THBS2","COL14A1","SFRP2","RBFOX1",
"ITGA11","ZNF385B","CHAD","AC007952.4","WWP2","CSGALNACT1","CTSC","NEGR1","NTNG1","DCN","KCTD8","CDH18","EGFL6","NTM","CCN3",
"COL6A6","THBS2","PTN","COL1A1","KHDRBS2","COL1A2","KCNIP4","COL1A1","IGF1","DKK2","CCDC80","THBS4","TNMD","FSTL5","FGF14",
"LINC02008","HBG2","EEF1A1","GYPA","KLF1","RPS8","RIPK4","BMP2","RIMS2","NRCAM","BMP6","SORBS2","FMOD","CYTL1","FRZB","CPE",
"PBK","CDK1","CDCA3","AURKB","UBE2C","RGS5","TRPC6","GUCY1A1","IGFBP7","AVPR1A","KLHL41","ACTC1","MYOG","RYR1","COBL","PECAM1","FLT1",
"ESAM","KDR","PLVAP","CRABP1","LINC01133","DCN","TMSB10","LINC01738","CXCL8","LAPTM5","KYNU","IL1B","PIK3R5","SLC35F1","MPZ","XKR4",
"CDH6","FOXD3")

DoHeatmap(subset(IL_57, downsample = 100), features = features, size = 3)
              



#Upstream targets of PTH1R             
#HAR1 IGF1           
CoveragePlot(IL_57, region = 'chr12-102454308-102459311', features = 'IGF1', assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE)
#HAR2 IGF1
CoveragePlot(IL_57, region = 'chr12-102470334-102472043', features = 'IGF1', assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE)

#IHH 
CoveragePlot(IL_57, region = 'chr2-219052800-219062546', features = 'IHH', assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE)

#SOX9
CoveragePlot(IL_57, region = 'chr17-72119670-72127765', features = 'SOX9', assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE)

#WNT7A
CoveragePlot(IL_57, region = 'chr3-13800305-13896025', features = 'SOX9', assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE)


#ZNF521 HAR1
chr18:25134105-25136383
CoveragePlot(IL_57, region = 'chr18-25077227-25078700', features = 'ZNF521', assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE)

#ZNF521 HCONDEL1
CoveragePlot(IL_57, region = 'chr18-25171076-25172079', features = 'ZNF521', assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE)

#ZNF521 HCONDEL2
CoveragePlot(IL_57, region = 'chr18-25175387-25176337', features = 'ZNF521', assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE)

#ZNF521 Hcondel3
CoveragePlot(IL_57, region = 'chr18-25369025-25369937', features = 'ZNF521', assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE)


#FOXP1 Locus

#FOXP1 HAR1
CoveragePlot(IL_57, region = 'chr3-70729363-70730582', features = 'FOXP1', assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE)

#FOXP1 HCondel 1
CoveragePlot(IL_57, region = 'chr3-71005490-71006659', features = 'FOXP1', assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE)
chr3:71005490-71006659

#FOXP1 Hcondel 2
CoveragePlot(IL_57, region = 'chr3-71055685-71057040', features = 'FOXP1', assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE)
chr3:71055685-71057040

#FOXP1 Hcondel 3
CoveragePlot(IL_57, region = 'chr3-71066641-71067561', features = 'FOXP1', assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE)
chr3:71066641-71067561

#FOXP1 HAR 2
CoveragePlot(IL_57, region = 'chr3-71303673-71305808', features = 'FOXP1', assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE)
chr3:71303673-71305808

#FOXP1 HAR 3
CoveragePlot(IL_57, region = 'chr3-71526081-71528694', features = 'FOXP1', assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE)
chr3:71526081-71528694


#FOXP2 HAR1
CoveragePlot(IL_57, region = 'chr7-114345191-114346456', features = 'FOXP2', assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE)

#FOXP2 HCONDEL1
CoveragePlot(IL_57, region = 'chr7-114696185-114710009', features = 'FOXP2', assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE)

#RUNX2- putative enhancer HAR1
CoveragePlot(IL_57, region = 'chr6-44916972-44919316', features = 'RUNX2', assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE)

#RUNX2 HCONDEL1
CoveragePlot(IL_57, region = 'chr6-45020513-45022590', features = 'RUNX2', assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE)

#RUNX2 HCONDEL2
CoveragePlot(IL_57, region = 'chr6-45771590-45773992', features = 'RUNX2', assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE)


chr17:72150834-72152468
CoveragePlot(IL_57, region = 'chr17-72150834-72152468', features = 'SOX9', assay = 'ATAC', expression.assay = 'SCT', peaks = TRUE)
