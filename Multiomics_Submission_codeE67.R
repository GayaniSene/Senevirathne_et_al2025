#this is the code for E67_Multiomics_IL



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
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(presto)
library("ggseqlogo")


#this is the code you would do before integrating spatial and single-cell data. To name the clusters etc. 

# the 10x hdf5 file contains both data types. 
inputdata.10x <- Read10X_h5("filtered_feature_bc_matrix.h5")

# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

# Create Seurat object
IL_67 <- CreateSeuratObject(counts = rna_counts)
IL_67[["percent.mt"]] <- PercentageFeatureSet(IL_67, pattern = "^MT-")

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

IL_67[["ATAC"]] <- chrom_assay


#We perform basic QC based on the number of detected molecules for each modality as well as mitochondrial percentage.

VlnPlot(IL_67, features = c("nCount_ATAC", "nCount_RNA", "percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()

#If you want to subset the dataset based on the QC values,
IL_67 <- subset(
  x = IL_67,
  subset = 
    percent.mt < 20
)


# RNA analysis
DefaultAssay(IL_67) <- "RNA"
IL_67 <- SCTransform(IL_67, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(IL_67) <- "ATAC"
IL_67 <- RunTFIDF(IL_67)
IL_67 <- FindTopFeatures(IL_67, min.cutoff = 'q0')
IL_67 <- RunSVD(IL_67)
IL_67 <- RunUMAP(IL_67, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

#We calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. We use this graph for UMAP visualization and clustering

IL_67 <- FindMultiModalNeighbors(IL_67, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
IL_67 <- RunUMAP(IL_67, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
IL_67 <- FindClusters(IL_67, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

#clustering based on gene expression, ATAC-seq, or WNN analysis.
p1 <- DimPlot(IL_67, reduction = "umap.rna", label = FALSE, label.size = 3.5, pt.size = 1.4, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(IL_67, reduction = "umap.atac",  label = FALSE, label.size = 3.5, pt.size = 1.4, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(IL_67, reduction = "wnn.umap",  label = TRUE, label.size = 3.5, pt.size = 2, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

p3


#Define the major groups of cell populations

DefaultAssay(IL_67) <- "SCT"
# we can also identify alternative protein and RNA markers for this cluster through
# differential expression

rna_markers0 <- FindMarkers(IL_67, ident.1 = 0, assay = "SCT")
head(rna_markers0)
head(rna_markers0, n = 50)


rna_markers1 <- FindMarkers(IL_67, ident.1 = 1, assay = "SCT")
head(rna_markers1)
head(rna_markers1, n = 50)

rna_markers2 <- FindMarkers(IL_67, ident.1 = 2, assay = "SCT")
head(rna_markers2)
head(rna_markers2, n = 50)

rna_markers3 <- FindMarkers(IL_67, ident.1 = 3, assay = "SCT")
head(rna_markers3)
head(rna_markers3, n = 50)

rna_markers4 <- FindMarkers(IL_67, ident.1 = 4, assay = "SCT")
head(rna_markers4)
head(rna_markers4, n = 50)

rna_markers5 <- FindMarkers(IL_67, ident.1 = 5, assay = "SCT")
head(rna_markers5)
head(rna_markers5, n = 50)

rna_markers6 <- FindMarkers(IL_67, ident.1 = 6, assay = "SCT")
head(rna_markers6)
head(rna_markers6, n = 50)

rna_markers7 <- FindMarkers(IL_67, ident.1 = 7, assay = "SCT")
head(rna_markers7)
head(rna_markers7, n = 50)

rna_markers8 <- FindMarkers(IL_67, ident.1 = 8, assay = "SCT")
head(rna_markers8)
head(rna_markers8, n = 50)

rna_markers9 <- FindMarkers(IL_67, ident.1 = 9, assay = "SCT")
head(rna_markers9)
head(rna_markers9, n = 50)

rna_markers10 <- FindMarkers(IL_67, ident.1 = 10, assay = "SCT")
head(rna_markers10)
head(rna_markers10, n = 50)

rna_markers11 <- FindMarkers(IL_67, ident.1 = 11, assay = "SCT")
head(rna_markers11)
head(rna_markers11, n = 50)

rna_markers12 <- FindMarkers(IL_67, ident.1 = 12, assay = "SCT")
head(rna_markers12)
head(rna_markers12, n = 50)

rna_markers13 <- FindMarkers(IL_67, ident.1 = 13, assay = "SCT")
head(rna_markers13)
head(rna_markers13, n = 50)

rna_markers14 <- FindMarkers(IL_67, ident.1 = 14, assay = "SCT")
head(rna_markers14)
head(rna_markers14, n = 50)

rna_markers15 <- FindMarkers(IL_67, ident.1 = 15, assay = "SCT")
head(rna_markers15)
head(rna_markers15, n = 50)

rna_markers16 <- FindMarkers(IL_67, ident.1 = 16, assay = "SCT")
head(rna_markers16)
head(rna_markers16, n = 50)



#Renaming the clusters

IL_67 <- RenameIdents(IL_67, '0' = 'RestingChondro','1' = 'Mes1','2' = 'Perichondrium')
IL_67 <- RenameIdents(IL_67, '3' = 'Proliferating1', '4' ='Mes2', '5' = 'Proliferating2')
IL_67 <- RenameIdents(IL_67, '6' = 'MyoProg+7','7' = 'Proliferating3', '8' = 'Proliferating4')
IL_67 <- RenameIdents(IL_67, '9' = 'TransMes', '10' = 'Osteoblasts', '11' = 'Perimysium')
IL_67 <- RenameIdents(IL_67, '12' = 'PrimaryErythro', '13' = 'Fibro1', '14' = 'Macrophages')
IL_67 <- RenameIdents(IL_67, '15' = 'Myocytes', '16' = 'ArterialEndo')
IL_67$celltype <- Idents(IL_67)


p1 <- DimPlot(IL_67, reduction = "umap.rna", label = FALSE, label.size = 3.5, pt.size = 1.4, repel = TRUE, cols = c('RestingChondro'='#055C9D', 'Mes1'= '#814827', 'Perichondrium' = '#60A3D9', 'Proliferating1'= '#0E86D4', 'Mes2'= '#A47551','Proliferating2'= '#0E86D4', 
                                                                                                                   'MyoProg+7'='#D43790', 'Proliferating3' =  '#0E86D4', 'Proliferating4' =  '#189AB4','TransMes'='#D2A87E', 
                                                                                                                   'Osteoblasts' =  '#BA0F30', 'Perimysium'= '#FF0080','PrimaryErythro'= '#F5631A',  'Fibro1' = '#2C5E1A', 'Macrophages'= '#603F8B','Myocytes'='#FB4570',
                                                                                                                   'ArterialEndo' = '#E55B13')) + ggtitle("RNA")
p2 <- DimPlot(IL_67, reduction = "umap.atac",  label = FALSE, label.size = 3.5, pt.size = 1.4, repel = TRUE, cols = c('RestingChondro'='#055C9D', 'Mes1'= '#814827', 'Perichondrium' = '#60A3D9', 'Proliferating1'= '#0E86D4', 'Mes2'= '#A47551','Proliferating2'= '#0E86D4', 
                                                                                                                     'MyoProg+7'='#D43790', 'Proliferating3' =  '#0E86D4', 'Proliferating4' =  '#189AB4','TransMes'='#D2A87E', 
                                                                                                                     'Osteoblasts' =  '#BA0F30', 'Perimysium'= '#FF0080','PrimaryErythro'= '#F5631A',  'Fibro1' = '#2C5E1A', 'Macrophages'= '#603F8B','Myocytes'='#FB4570',
                                                                                                                     'ArterialEndo' = '#E55B13')) + ggtitle("ATAC")
p3 <- DimPlot(IL_67, reduction = "wnn.umap",  label = TRUE, label.size = 3.5, pt.size = 2, repel = TRUE, cols = c('RestingChondro'='#055C9D', 'Mes1'= '#814827', 'Perichondrium' = '#60A3D9', 'Proliferating1'= '#0E86D4', 'Mes2'= '#A47551','Proliferating2'= '#0E86D4', 
                                                                                                                 'MyoProg+7'='#D43790', 'Proliferating3' =  '#0E86D4', 'Proliferating4' =  '#189AB4','TransMes'='#D2A87E', 
                                                                                                                 'Osteoblasts' =  '#BA0F30', 'Perimysium'= '#FF0080','PrimaryErythro'= '#F5631A',  'Fibro1' = '#2C5E1A', 'Macrophages'= '#603F8B','Myocytes'='#FB4570',
                                                                                                                 'ArterialEndo' = '#E55B13')) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))


saveRDS(IL_67, file = "IL_673k_final.rds")

#marker gene plots in Extended Data Fig3. 
DefaultAssay(IL_67) <- "SCT"
features <- c("MATN4","CNMD","UCMA","COL9A2","MATN3","COMP","SNORC","SHOX","DLX5", "RUNX2", "FGFBP2","IHH","THBS2","COL2A1","COL1A2","COL10A1", "PAX3", "PAX7", "TNMD", "IGFBP7", "MPZ", "SP7", "SOX9", "VEGFA", "VEGFB", "ZNF521", "PTH1R","FOXD3")

DotPlot(IL_67, features = features) + RotatedAxis()


features <- c("LSAMP","UCMA","MT-CO1","GLIS3","KCNMA1", "SCN7A","NTNG1","COL6A6","ADAMTSL1","EPHA7",
              "THBS2","SFRP2","COL1A1","COL1A2","COL12A1","PKIA","GLIS3","THRB","ANKH","ANKRD28",
              "MGAT4C","DCLK1","NTRK2","COL1A1","CXCL12","COL2A1","SOX6","AC007952.4","COL9A1","ACAN",
              "VGLL3","PDLIM3","PDE1C","PAX7","MYF5","TNC","ERG","PLCB1","CYTL1","MALAT1","COL2A1","CNMD","MATN1","COL9A2",
              "TOP2A","DIAPH3","CENPF","MKI67","HIST1H1A","NR4A2","NR4A3","BMP2","CHMP1B","BMP6","TNMD","FSTL5","FGF14","LINC02008","THBS4",
              "GYPB","RHAG","SPTA1","GYPA","FNDC3B","CXCL8","CCL3","KYNU","IL1B","LCP2","KLHL41","ARPP21","SYN2","MYOG","TMCC3",
              "ADGRF5","MYCT1","CD93","CDH5","DIPK2B")
              
                            

DoHeatmap(subset(IL_67, downsample = 100), features = features, size = 3)







