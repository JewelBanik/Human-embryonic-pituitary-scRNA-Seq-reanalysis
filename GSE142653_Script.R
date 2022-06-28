#GSE142653

metadata.File <- read.csv("GSE142653_pit_dev_CellInfo.csv", header = T, row.names = 1)

counts.File <- read.csv("GSE142653_pit_dev_5181_count.csv", header = T, row.names = 1)

mm <- match(rownames(metadata.File), colnames(counts.File))

matched.File <- counts.File[mm]

library(Seurat)
library(SeuratObject)
library(patchwork)
library(ggplot2)
library(dplyr)

data <- CreateSeuratObject(counts = matched.File, meta.data = metadata.File)

VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
#saveRDS(data, file = "GSE142653.rds")

plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1+plot2

GSE142653 <- NormalizeData(data, normalization.method = "LogNormalize", 
                           scale.factor = 100000)
GSE142653 <- FindVariableFeatures(data, selection.method = "vst", 
                                  nfeatures = 2000)

all.genes <- rownames(GSE142653)
GSE142653 <- ScaleData(GSE142653, features = all.genes)

GSE142653 <- RunPCA(GSE142653, features = VariableFeatures(object = GSE142653))
DimPlot(GSE142653, reduction = "pca")

GSE142653 <- JackStraw(GSE142653, num.replicate = 100)
GSE142653 <- ScoreJackStraw(GSE142653, dims = 1:20)
JackStrawPlot(GSE142653, dims = 1:20)
ElbowPlot(GSE142653)


GSE142653 <- FindNeighbors(GSE142653, dims = 1:20)
GSE142653 <- FindClusters(GSE142653, resolution = 1)

GSE142653 <- RunUMAP(GSE142653, dims = 1:20)
DimPlot(GSE142653, reduction = 'umap')



