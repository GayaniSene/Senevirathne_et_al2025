#This is the code for E72_IL_multiomics

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
IL_72 <- CreateSeuratObject(counts = rna_counts)
IL_72[["percent.mt"]] <- PercentageFeatureSet(IL_72, pattern = "^MT-")

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

IL_72[["ATAC"]] <- chrom_assay


#We perform basic QC based on the number of detected molecules for each modality as well as mitochondrial percentage.

VlnPlot(IL_72, features = c("nCount_ATAC", "nCount_RNA", "percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()

#If you want to subset the dataset based on the QC values,
IL_72 <- subset(
  x = IL_72,
  subset = 
    percent.mt < 10
)


# RNA analysis
DefaultAssay(IL_72) <- "RNA"
IL_72 <- SCTransform(IL_72, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(IL_72) <- "ATAC"
IL_72 <- RunTFIDF(IL_72)
IL_72 <- FindTopFeatures(IL_72, min.cutoff = 'q0')
IL_72 <- RunSVD(IL_72)
IL_72 <- RunUMAP(IL_72, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

#We calculate a WNN graph, representing a weighted combination of RNA and ATAC-seq modalities. We use this graph for UMAP visualization and clustering

IL_72 <- FindMultiModalNeighbors(IL_72, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
IL_72 <- RunUMAP(IL_72, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
IL_72 <- FindClusters(IL_72, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

#clustering based on gene expression, ATAC-seq, or WNN analysis.
p1 <- DimPlot(IL_72, reduction = "umap.rna", label = FALSE, label.size = 3.5, pt.size = 1.4, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(IL_72, reduction = "umap.atac",  label = FALSE, label.size = 3.5, pt.size = 1.4, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(IL_72, reduction = "wnn.umap",  label = TRUE, label.size = 3.5, pt.size = 2, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

p3


#Define the major groups of cell populations

DefaultAssay(IL_72) <- "SCT"
# we can also identify alternative protein and RNA markers for this cluster through
# differential expression

rna_markers0 <- FindMarkers(IL_72, ident.1 = 0, assay = "SCT")
head(rna_markers0)
head(rna_markers0, n = 50)


rna_markers1 <- FindMarkers(IL_72, ident.1 = 1, assay = "SCT")
head(rna_markers1)
head(rna_markers1, n = 50)

rna_markers2 <- FindMarkers(IL_72, ident.1 = 2, assay = "SCT")
head(rna_markers2)
head(rna_markers2, n = 50)

rna_markers3 <- FindMarkers(IL_72, ident.1 = 3, assay = "SCT")
head(rna_markers3)
head(rna_markers3, n = 50)

rna_markers4 <- FindMarkers(IL_72, ident.1 = 4, assay = "SCT")
head(rna_markers4)
head(rna_markers4, n = 50)

rna_markers5 <- FindMarkers(IL_72, ident.1 = 5, assay = "SCT")
head(rna_markers5)
head(rna_markers5, n = 50)

rna_markers6 <- FindMarkers(IL_72, ident.1 = 6, assay = "SCT")
head(rna_markers6)
head(rna_markers6, n = 50)

rna_markers7 <- FindMarkers(IL_72, ident.1 = 7, assay = "SCT")
head(rna_markers7)
head(rna_markers7, n = 50)

rna_markers8 <- FindMarkers(IL_72, ident.1 = 8, assay = "SCT")
head(rna_markers8)
head(rna_markers8, n = 50)

rna_markers9 <- FindMarkers(IL_72, ident.1 = 9, assay = "SCT")
head(rna_markers9)
head(rna_markers9, n = 50)

rna_markers10 <- FindMarkers(IL_72, ident.1 = 10, assay = "SCT")
head(rna_markers10)
head(rna_markers10, n = 50)

rna_markers11 <- FindMarkers(IL_72, ident.1 = 11, assay = "SCT")
head(rna_markers11)
head(rna_markers11, n = 50)

rna_markers12 <- FindMarkers(IL_72, ident.1 = 12, assay = "SCT")
head(rna_markers12)
head(rna_markers12, n = 50)

rna_markers13 <- FindMarkers(IL_72, ident.1 = 13, assay = "SCT")
head(rna_markers13)
head(rna_markers13, n = 50)

rna_markers14 <- FindMarkers(IL_72, ident.1 = 14, assay = "SCT")
head(rna_markers14)
head(rna_markers14, n = 50)

rna_markers15 <- FindMarkers(IL_72, ident.1 = 15, assay = "SCT")
head(rna_markers15)
head(rna_markers15, n = 50)

rna_markers16 <- FindMarkers(IL_72, ident.1 = 16, assay = "SCT")
head(rna_markers16)
head(rna_markers16, n = 50)

rna_markers17 <- FindMarkers(IL_72, ident.1 = 17, assay = "SCT")
head(rna_markers17)
head(rna_markers17, n = 50)

rna_markers18 <- FindMarkers(IL_72, ident.1 = 18, assay = "SCT")
head(rna_markers18)
head(rna_markers18, n = 50)

rna_markers19 <- FindMarkers(IL_72, ident.1 = 19, assay = "SCT")
head(rna_markers19)
head(rna_markers19, n = 50)

rna_markers20 <- FindMarkers(IL_72, ident.1 = 20, assay = "SCT")
head(rna_markers20)
head(rna_markers20, n = 50)




#Renaming the clusters

IL_72 <- RenameIdents(IL_72, '0' = 'RestingChondro','1' = 'Mes1','2' = 'Myocytes')
IL_72 <- RenameIdents(IL_72, '3' = 'Proliferating1', '4' ='Perichondrium', '5' = 'MyoProg')
IL_72 <- RenameIdents(IL_72, '6' = 'Perimysium','7' = 'TransMes', '8' = 'Mes2')
IL_72 <- RenameIdents(IL_72, '9' = 'Proliferating2', '10' = 'Fibro1', '11' = 'TransMes')
IL_72 <- RenameIdents(IL_72, '12' = 'Myocytes', '13' = 'Perimysium2', '14' = 'PrimaryErythro')
IL_72 <- RenameIdents(IL_72, '15' = 'Erythrocytes', '16' = 'VenousEndo', '17' = 'Macrophages')
IL_72 <- RenameIdents(IL_72, '18' = 'Myoprog+Pax7', '19'= 'Schwann', '20'='Neurons')
IL_72$celltype <- Idents(IL_72)




p1 <- DimPlot(IL_72, reduction = "umap.rna", label = FALSE, label.size = 3.5, pt.size = 1.4, repel = TRUE, cols = c('RestingChondro'='#055C9D', 'Mes1'= '#814827', 'Myocytes'='#FB4570', 'Proliferating1'= '#0E86D4', 'Perichondrium' = '#60A3D9', 'MyoProg'='#D43790',
                                                                                                                   'Perimysium'= '#FF0080', 'TransMes'= '#D2A87E' , 'Mes2'= '#A47551', 'Proliferating2'= '#0E86D4', 'Fibro1' = '#2C5E1A', 'TransMes'='#D2A87E',
                                                                                                                   'Myocytes' = '#FB4570', 'Perimysium2'= '#FF0080', 'PrimaryErythro'= '#F5631A', 'Erythrocytes'= '#FD7F20', 'VenousEndo'= '#D85E00',
                                                                                                                   'Macrophages'= '#603F8B','Myoprog+Pax7'='#D43790','Schwann'='#FFDD64','Neurons'='#FFDD64')) + ggtitle("RNA")
p2 <- DimPlot(IL_72, reduction = "umap.atac",  label = FALSE, label.size = 3.5, pt.size = 1.4, repel = TRUE, cols = c('RestingChondro'='#055C9D', 'Mes1'= '#814827', 'Myocytes'='#FB4570', 'Proliferating1'= '#0E86D4', 'Perichondrium' = '#60A3D9', 'MyoProg'='#D43790',
                                                                                                                     'Perimysium'= '#FF0080', 'TransMes'= '#D2A87E' , 'Mes2'= '#A47551', 'Proliferating2'= '#0E86D4', 'Fibro1' = '#2C5E1A', 'TransMes'='#D2A87E',
                                                                                                                     'Myocytes' = '#FB4570', 'Perimysium2'= '#FF0080', 'PrimaryErythro'= '#F5631A', 'Erythrocytes'= '#FD7F20', 'VenousEndo'= '#D85E00',
                                                                                                                     'Macrophages'= '#603F8B','Myoprog+Pax7'='#D43790','Schwann'='#FFDD64','Neurons'='#FFDD64')) + ggtitle("ATAC")
p3 <- DimPlot(IL_72, reduction = "wnn.umap",  label = TRUE, label.size = 3.5, pt.size = 2, repel = TRUE, cols = c('RestingChondro'='#055C9D', 'Mes1'= '#814827', 'Myocytes'='#FB4570', 'Proliferating1'= '#0E86D4', 'Perichondrium' = '#60A3D9', 'MyoProg'='#D43790',
                                                                                                                 'Perimysium'= '#FF0080', 'TransMes'= '#D2A87E' , 'Mes2'= '#A47551', 'Proliferating2'= '#0E86D4', 'Fibro1' = '#2C5E1A', 'TransMes'='#D2A87E',
                                                                                                                 'Myocytes' = '#FB4570', 'Perimysium2'= '#FF0080', 'PrimaryErythro'= '#F5631A', 'Erythrocytes'= '#FD7F20', 'VenousEndo'= '#D85E00',
                                                                                                                 'Macrophages'= '#603F8B','Myoprog+Pax7'='#D43790','Schwann'='#FFDD64','Neurons'='#FFDD64')) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

saveRDS(IL_72, file = "IL_723k_final.rds")


p3 <- DimPlot(IL_72, reduction = "wnn.umap",  label = FALSE, label.size = 3.5, pt.size = 2, repel = TRUE, cols = c('RestingChondro'='#68598A', 'Mes1'= '#8D5876', 'Myocytes'='#B9BFC9', 'Proliferating1'= '#B9BFC9', 'Perichondrium' = '#B9BFC9', 'MyoProg'='#B9BFC9',
                                                                                                                 'Perimysium'= '#B9BFC9', 'TransMes'= '#B9BFC9' , 'Mes2'= '#B9BFC9', 'Proliferating2'= '#B9BFC9', 'Fibro1' = '#B9BFC9', 'TransMes'='#B9BFC9',
                                                                                                                 'Myocytes' = '#7280A9', 'Perimysium2'= '#B9BFC9', 'PrimaryErythro'= '#B9BFC9', 'Erythrocytes'= '#B9BFC9', 'VenousEndo'= '#B9BFC9',
                                                                                                                 'Macrophages'= '#603F8B','Myoprog+Pax7'='#B9BFC9','Schwann'='#B9BFC9','Neurons'='#B9BFC9')) + ggtitle("WNN")

#Import seurat clusters for the velocity analysis using scvelo
#Information from wnn.UMAP for indivdual cells are saved to .csv files. 

IL_72$barcode <- colnames(IL_72)
IL_72$WNN.UMAP_1 <- IL_72@reductions$wnn.umap@cell.embeddings[,1]
IL_72$WNN.UMAP_2 <- IL_72@reductions$wnn.umap@cell.embeddings[,2]
write.csv(IL_72@meta.data, file='metadata.csv', quote=F, row.names=F)

library(Matrix)
counts_matrix <- GetAssayData(IL_72, assay='SCT', slot='counts')

writeMM(counts_matrix, file=paste0('counts.mtx'))
write.csv(IL_72@reductions$pca@cell.embeddings, file='pca.csv', quote=F, row.names=F)

#cell_data table with annotations and barcodes for the input file for SCENIC+ 
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',
  quote=F,row.names=F,col.names=F
)

#to export barcodes and celltypes
library(tibble)
export_df <- IL_72@meta.data %>% 
  rownames_to_column("barcodes")

head(export_df)
write.csv(export_df, "cell_data.csv")



