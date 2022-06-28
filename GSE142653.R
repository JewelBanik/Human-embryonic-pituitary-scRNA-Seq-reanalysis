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
library(tidyverse)

# Here we are filtering out genes that are are expressed in 3 or fewer cells and are filtering cells with complexity less than 200 genes.
data <- CreateSeuratObject(counts = matched.File, meta.data = metadata.File, project = "projects")
data
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

gse <- NormalizeData(object = data, normalization.method = "LogNormalize", scale.factor = 100000)

gse <- FindVariableFeatures(gse, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(gse)

all.genes <- rownames(gse)
gse <- ScaleData(gse, features = all.genes)

#gse2 <- ScaleData(gse, features = all.genes, vars.to.regress = 'percent.mito')

gse <- RunPCA(gse, features = VariableFeatures(object = gse))
DimPlot(gse, reduction = 'pca', group.by = 'cell_type', label = T)

gse <- JackStraw(gse, num.replicate = 100)
gse <- ScoreJackStraw(gse, dims = 1:20)
JackStrawPlot(gse, dims = 1:20)

ElbowPlot(gse)


gse <- FindNeighbors(gse, dims = 1:20)
gse <- FindClusters(gse, resolution = 1)

gse <- RunUMAP(gse, dims = 1:20)

DimPlot(gse, reduction = 'umap', group.by = 'cell_type', label = T, pt.size = 1.5, repel = F)+
  labs(title = "Fetal Human Pituitary_GSE142653")

png(file = "UMAP_cell_types.png", height = 4, width = 8, res = 1200, units = 'in')
DimPlot(gse, reduction = 'umap', group.by = 'cell_type', label = T, pt.size = 1.5, 
        repel = F, label.size = 2) + labs(title = "Fetal Human Pituitary_Cell Types")
dev.off()


DimPlot(gse, reduction = 'umap', group.by = 'orig.ident', pt.size = 1.5)+
  labs(title = "Fetal Human Pituitary_GSE142653")

png(file = "UMAP_orig.ident.png", height = 4, width = 8, res = 1200, units = 'in')
DimPlot(gse, reduction = 'umap', group.by = 'orig.ident', label = F, pt.size = 1.5, 
        repel = F, label.size = 2, order = T) + labs(title = "Fetal Human Pituitary_Temporal Trajectory")
dev.off()

####DimPlot with color palette: 
colPallete <- c('gray0', 'gray10', 'gray20', 'gray30', 'gray40', 'gray50', 
                'deepskyblue', 'deepskyblue1', 'deepskyblue2', 'deepskyblue3',
                'dodgerblue', 'dodgerblue1', 'dodgerblue2', 'dodgerblue3', 
                'dodgerblue4', 'greenyellow', 'lawngreen', 'green', 'green1',
                'green2', 'green3')
# https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf

png(file = "UMAP_orig.ident_colPalette.png", height = 4, width = 8, res = 1200, units = 'in')
DimPlot(gse, reduction = 'umap', group.by = 'orig.ident', label = F, pt.size = 1.5, 
        repel = F, label.size = 2, order = T, cols = colPallete)+
  labs(title = "Fetal Human Pituitary_Temporal Trajectory")
dev.off()





DimPlot(gse, reduction = 'umap', group.by = 'subcluster', label = T, pt.size = 1, 
        repel = T, label.size = 3)+labs(title = "Fetal Human Pituitary_GSE142653")

png(file = "UMAP_subclusters.png", height = 4, width = 8, res = 1200, units = 'in')

DimPlot(gse, reduction = 'umap', group.by = 'subcluster', label = T, pt.size = 1.5, 
        repel = T, label.size = 2, cells = exAPcells) + labs(title = "Fetal Human Pituitary_Cell Subclusters")+
  guides(guides(color = guide_legend(override.aes = list(size=1), ncol=2)))

dev.off()


features <- c("PITX1", "PITX2", "SOX2", "PROP1", "POU1F1", "TBX19", "NR5A1", "POMC", 
              'GH1', 'GH2', 'TSHB', 'PRL', 'FSHB', 'LHB', 'PDGFRA', 'OTX2', 'LHX2', 
              'RAX', 'COL25A1', 'PECAM1', 'HBQ1', 'PTPRC', 'MSI1', 'MSI2')


n_cells <- FetchData(gse, vars = c('ident', 'cell_type', 'orig.ident'))
n <- n_cells %>%
dplyr::count(cell_type, orig.ident) %>%
tidyr::spread(orig.ident, n, fill = 0)

View(n)



rs <- rowSums(n[ ,2:22])
cell.sum <- cbind(n, Sum=rs)
View(cell.sum)

cs <- colSums(cell.sum[2:23])
cs <- c("Sum", cs)#to add colsums to the dataframe as rbind. 

cell.sum <- rbind(cell.sum, cs)
View(cell.sum)

library(stargazer)
stargazer(cell.sum, type = "text", summary = F, out = "CellSum.txt")

DotPlot(gse, features = c("TBX19", 'PAX7', 'POMC', 'GH2', 'NR5A1', 'FSHB', 'LHB', 'GATA2', 'TSHB', 'PRL'), group.by = "orig.ident")
DotPlot(gse, features = c("TBX19",'POMC', 'GH2', 'NR5A1', 'FSHB', 'LHB', 'GATA2', 'TSHB', 'PRL'), group.by = "cell_type")

png(file = "DotPlot_canonical endo markers_cell_types.png", height = 4, width = 10, res = 1200, units = 'in')
DotPlot(gse, features = c("TBX19",'POMC', 'GH2', 'NR5A1', 'FSHB', 'LHB', 'GATA2', 'TSHB', 'PRL'), group.by = "cell_type")
dev.off()

png(file = "DotPlot_canonical lineage markers_cell_types.png", height = 4, width = 10, res = 1200, units = 'in')
DotPlot(gse, features = c("PITX1", "PROP1", "POU1F1", "SOX2",'OTX2', 'PDGFRA', 'HBQ1', 'PECAM1'), group.by = "cell_type")
dev.off()


VlnPlot(gse, features = 'GAPDH', group.by = 'orig.ident')

FeaturePlot(gse2, features = c("PITX1", "SOX2", "PROP1", "POU1F1", "TBX19", "NR5A1", "POMC",'GH2', 
                               'TSHB', 'PRL', 'FSHB', 'LHB', 'PDGFRA', 'OTX2','RAX','PECAM1', 'HBQ1', 
                               'PTPRC', 'MSI1', 'MSI2'), cols = c('blue', 'red'), ncol = 5)

markers <- FindAllMarkers(gse, min.pct = 0.25, logfc.threshold = 0.25)#go with this one as standard since 
#most of the canonical markers for each cell type are represented by this min.pct threshold. 

top5 <- markers %>%
group_by(cluster)%>%
top_n(n=5, wt=avg_log2FC)

View(top5)

markers2 <- FindAllMarkers(gse2, min.pct = 0.1, logfc.threshold = 0.25)

#This plot shows that folliculostellate cells (FSC) were not present in human pit at P25, 
#marked by no cells that express these 4 genes. 
VlnPlot(gse, features = c("S100B", 'FXYD1', 'GSTM2', 'CAPN6'), group.by = 'cell_type', ncol = 1)
#This plot shows the expression of PROP1 regulators in the Prop1 expressing cells
VlnPlot(gse, features = c("ZFP36L1", 'ANXA1', 'NFIB', 'ZNF521', 'NR2F2'), group.by = 'cell_type', ncol = 1, pt.size =0.2)

###This part of code explains how I isolated only Endocrine cells from the whole pit data

endo<- c("Stem", "Corticotrope", 'Pro.PIT1', 'CC', 'Thyrotrope', 'Pre.Gonado', 'Somatotrope', 'Lactotrope', 'Gonadotrope')

result <- filter(gse@meta.data, grepl(paste(endo, collapse="|"), cell_type)) #this line returns only 
#endocrine cells metadata; however, this vector can't be embedded in DimPlot. 

APcells<- grep(paste(endo, collapse = "|"), gse@meta.data$cell_type, value = F) #Value = if FALSE, 
#a vector containing the (integer) indices of the matches determined by grep is returned, 
#and if TRUE, a vector containing the matching elements themselves is returned.I needed only the integer
#indices, thus value = F. Here, I'm isolating only endocrine cells (n=2388). 

DimPlot(object = gse, reduction = 'umap', cells = APcells, group.by = 'cell_type', pt.size = 1.5)
DimPlot(object = gse, reduction = 'umap', cells = APcells, group.by = 'orig.ident', pt.size = 1.5)
#Now you see, here some of the endocrine cells (20) that are scattered throught the clusters and seems like
#they don't represent a distinct cell type. Even though their PC_1 values are higher than other cells,
#I don't think they influence the canonical clusters significantly; hence, I preferred to exclude them 
#from the plot, at least for the sake of aesthetic of the plot. 


plot<- DimPlot(object = gse, reduction = 'umap', cells = APcells, group.by = 'orig.ident', pt.size = 1.5)
HoverLocator(plot = plot, information = FetchData(gse, vars = c("ident", "PC_1", "nFeature_RNA", "APcells"))) #Here, 
#in this hoverLocator, I have manually noted the indices of outliers and excluded them from the plot. 
ex<- grep(pattern = paste(c(2246, 3318, 4885, 2997, 2414, 4318, 5201, 7635, 2522, 6636, 6996, 4888, 5024,
                            3552, 6155, 2572, 2738, 3153, 3627, 4320), collapse = '|'), 
          x = gse@meta.data$nFeature_RNA, value = F)
# Here, first, I searched cells according to nFeature_RNA values got from the hoverLocator of the cells 
#that I want to exclude. It returned 33 features. But, I provided 20 inputs. So, I have to figure out 
#the cells with the same number of features. 

rnames<- row.names.data.frame(gse@meta.data[ex, ]) # I stored the cells according to returned cell indices
View(cbind(rnames, ex)) # here, I matched the cell names from the hoverLocator and noted their indices (20).

#Now delete these 20 cells based on their indices from APcells. 
exIndex<- c(2, 15, 830, 1307, 3246, 3241, 3257, 2853, 4028, 3531, 2331, 2342, 2353, 3506, 1353, 2304, 
             2306, 3427, 2678, 2423)

exAPcells<- APcells[! APcells %in% exIndex]
length(exAPcells) # should return n=2368 cells. 20 outlier cells are deleted. 

#Now, create a DimPlot only with Endocrine cells. YahooOO!!!
DimPlot(object = gse, reduction = 'umap', cells = exAPcells, group.by = 'orig.ident', pt.size = 1.5)

png(file = "UMAP_endo_orig.ident.png", height = 4, width = 6.5, res = 1200, units = 'in')
DimPlot(gse, reduction = 'umap', group.by = 'orig.ident', cells = exAPcells, label = F, pt.size = 1.5, 
        repel = F, label.size = 2) + labs(title = "Fetal Human Pituitary_Endocrine Lineage Cells_Timepoint")
dev.off()



DimPlot(object = gse, reduction = 'umap', cells = exAPcells, group.by = 'cell_type', pt.size = 1.5)+
  labs(title = "Endocrine lineage cells of fetal human pituitary")

png(file = "UMAP_endo_cell_types.png", height = 4, width = 6, res = 1200, units = 'in')
DimPlot(gse, reduction = 'umap', group.by = 'cell_type', cells = exAPcells, label = T, pt.size = 1.5, 
        repel = F, label.size = 2) + labs(title = "Fetal Human Pituitary_Endocrine Lineage Cells_Cell Types")
dev.off()

  #coord_flip()+
  #coord_equal()+
  #coord_trans()+
  #coord_polar()+
  #coord_fixed()


DimPlot(object = gse, reduction = 'umap', cells = exAPcells, group.by = 'cell_type', pt.size = 1.5)+coord_flip()
VlnPlot(object = gse, features = c("LSM14A", 'LSM14B'), group.by = 'cell_type', same.y.lims = T, ncol = 1)
VlnPlot(object = gse, features = c('MSI1', 'MSI2'), group.by = 'cell_type', same.y.lims = T, ncol = 1)#this 
#plot implies that the expression level of MSI1 and MSI2 is kind of the same (y-axis) but it seems like 
#more cells express MSI1 than MSI2. 
VlnPlot(object = gse, features = c("LSM14A", 'LSM14B'), group.by = 'orig.ident', same.y.lims = T, ncol = 1)
VlnPlot(object = gse, features = c('MSI1', 'MSI2'), group.by = 'orig.ident', same.y.lims = T, ncol = 1)


tiff(file = "MSI1_MSI2_vlnPlot_orig.ident.tiff", height = 4, width = 6, res = 300, units = 'in', pointsize = 6,compression = 'lzw')

png(file = "vlnPlot_MSI1_MSI2_orig.ident.png", height = 6, width = 6, res = 1200, units = 'in')
VlnPlot(object = gse, features = c('MSI1', 'MSI2'), group.by = 'orig.ident', same.y.lims = T, ncol = 1)
dev.off()


# jpeg(file = "MSI1_MSI2_vlnPlot_orig.ident.jpg", height = 6, width = 6, res = 1200, units = 'in')
# VlnPlot(object = gse, features = c('MSI1', 'MSI2'), group.by = c('cell_type'), same.y.lims = T, ncol = 1)
# dev.off()

png(file = "MSI1_MSI2_vlnPlot_cell_types.png", height = 6, width = 6, res = 1200, units = 'in')
VlnPlot(object = gse, features = c('MSI1', 'MSI2'), group.by = 'cell_type', same.y.lims = T, ncol = 1)
dev.off()

png(file = "vlnPlot_LSM14A_LSM14B_cell_types.png", height = 6, width = 6, res = 1200, units = 'in')
VlnPlot(object = gse, features = c('LSM14A', 'LSM14B'), group.by = 'cell_type', same.y.lims = T, ncol = 1)
dev.off()

png(file = "vlnPlot_LSM14A_LSM14B_orig.ident.png", height = 6, width = 6, res = 1200, units = 'in')
VlnPlot(object = gse, features = c('LSM14A', 'LSM14B'), group.by = 'orig.ident', same.y.lims = T, ncol = 1)
dev.off()





#library(Cairo)
# Cairo::Cairo(
#   30, #length
#   30, #width
#   file = paste("nameofplot", ".png", sep = ""),
#   type = "png", #tiff
#   bg = "transparent", #white or transparent depending on your requirement 
#   dpi = 300,
#   units = "cm" #you can change to pixels etc 
# )
# plot(p) #p is your graph object 
# dev.off()


VlnPlot(object = gse2, features = c("LSM14A", 'LSM14B'), group.by = 'cell_type', same.y.lims = T)
FeaturePlot(gse, features = c('LSM14A', "LSM14B"), cols = c('yellow', 'red'), pt.size = 1.5, cells = exAPcells, order = T, min.cutoff = 'q10')
FeaturePlot(gse, features = c('MSI1', "MSI2"), cols = c('grey', 'blue'), pt.size = 1.5, cells = exAPcells, min.cutoff = 'q10', order = T)

png(file = "FeaturePlot_endo_MSI1_MSI2.png", height = 3.5, width = 8, res = 1200, units = 'in')
FeaturePlot(gse, features = c('MSI1', "MSI2"), cols = c('grey', 'blue'), pt.size = 1, 
            cells = exAPcells, min.cutoff = 'q10', order = T)
dev.off()


png(file = "FeaturePlot_MSI1_MSI2.png", height = 3.5, width = 8, res = 1200, units = 'in')
FeaturePlot(gse, features = c('MSI1', "MSI2"), cols = c('grey', 'blue'), pt.size = 1, min.cutoff = 'q10', order = T, coord.fixed = T)
dev.off()

png(file = "FeaturePlot_MSI1.png", height = 4, width = 8, res = 1200, units = 'in')
FeaturePlot(gse, features = 'MSI1', cols = c('grey', 'blue'), pt.size = 1, min.cutoff = 'q10', order = T, coord.fixed = T)
dev.off()

png(file = "FeaturePlot_MSI2.png", height = 4, width = 8, res = 1200, units = 'in')
FeaturePlot(gse, features = 'MSI2', cols = c('grey', 'blue'), pt.size = 1, min.cutoff = 'q10', order = T, coord.fixed = T)
dev.off()

png(file = "CombinePlot_MSI1_cell_type.png", height = 8.5, width = 10, res = 1200, units = 'in')
FeaturePlot(gse, features = 'MSI1', cols = c('grey', 'blue'), pt.size = 1.5, min.cutoff = 'q10', order = T, combine = T)+
  VlnPlot(object = gse, features = 'MSI1', group.by = 'cell_type', ncol = 1)
dev.off()


png(file = "CombinePlot_MSI2_cell_type.png", height = 8.5, width = 10, res = 1200, units = 'in')
FeaturePlot(gse, features = 'MSI2', cols = c('grey', 'blue'), pt.size = 1.5, min.cutoff = 'q10', order = T, combine = T)+
  VlnPlot(object = gse, features = 'MSI2', group.by = 'cell_type', ncol = 1)
dev.off()


png(file = "CombinePlot_MSI1_orig.ident.png", height = 8.5, width = 10, res = 1200, units = 'in')
FeaturePlot(gse, features = 'MSI1', cols = c('grey', 'blue'), pt.size = 1.5, min.cutoff = 'q10', order = T, combine = T)+
  VlnPlot(object = gse, features = 'MSI1', group.by = 'orig.ident', ncol = 1)
dev.off()


png(file = "CombinePlot_MSI2_orig.ident.png", height = 8.5, width = 10, res = 1200, units = 'in')
FeaturePlot(gse, features = 'MSI2', cols = c('grey', 'blue'), pt.size = 1.5, min.cutoff = 'q10', order = T, combine = T)+
  VlnPlot(object = gse, features = 'MSI2', group.by = 'orig.ident', ncol = 1)
dev.off()

png(file = "FeaturePlot_endo_MSI1.png", height = 4.5, width = 6, res = 1200, units = 'in')
FeaturePlot(gse, features = 'MSI1', cols = c('grey', 'blue'), pt.size = 1.5, 
            cells = exAPcells, min.cutoff = 'q10', order = T, coord.fixed = F)
dev.off()

png(file = "FeaturePlot_endo_MSI2.png", height = 4.5, width = 6, res = 1200, units = 'in')
FeaturePlot(gse, features = 'MSI2', cols = c('grey', 'blue'), pt.size = 1.5, 
            cells = exAPcells, min.cutoff = 'q10', order = T, coord.fixed = F)
dev.off()



png(file = "MSI1_MSI2_vlnPlot_cell_types.png", height = 6, width = 6, res = 1200, units = 'in')
VlnPlot(object = gse, features = c('MSI1', 'MSI2'), group.by = 'dd', same.y.lims = T, ncol = 1, )
dev.off()


dd<- FetchData(object = gse, vars= colnames(x = gse)[exAPcells])

ss<- subset(x = gse, subset = colnames(gse)[exAPcells])





plot1 <- DimPlot(gse, group.by = 'orig.ident')
plot2 <- FeatureScatter(gse, feature1 = "MSI1", feature2 = "MSI2", group.by = 'orig.ident')
# Combine two plots
plot1 + plot2



FeaturePlot(gse, features = c('MSI1', "MSI2"), cols = c('yellow', 'red'), pt.size = 1.5, cells = 1:170)
FeaturePlot(gse, features = c('MSI1', "MSI2"), cols = c('yellow', 'red'), pt.size = 1.5, cells = 171:347)
FeaturePlot(gse, features = c('MSI1', "MSI2"), cols = c('yellow', 'red'), pt.size = 1.5, cells = 348:561)
FeaturePlot(gse, features = c('MSI1', "MSI2"), cols = c('yellow', 'red'), pt.size = 1.5, cells = 562:695)
FeaturePlot(gse, features = c('MSI1', "MSI2"), cols = c('yellow', 'red'), pt.size = 1.5, cells = 696:866)
FeaturePlot(gse, features = c('MSI1', "MSI2"), cols = c('yellow', 'red'), pt.size = 1.5, cells = 867:995)
FeaturePlot(gse, features = c('MSI1', "MSI2"), cols = c('yellow', 'red'), pt.size = 1.5, cells = 996:1231)
FeaturePlot(gse, features = c('MSI1', "MSI2"), cols = c('yellow', 'red'), pt.size = 1.5, cells = 1232:1384)
FeaturePlot(gse, features = c('MSI1', "MSI2"), cols = c('yellow', 'red'), pt.size = 1.5, cells = 1385:1563)
FeaturePlot(gse, features = c('MSI1', "MSI2"), cols = c('yellow', 'red'), pt.size = 1.5, cells = 1564:1638)
FeaturePlot(gse, features = c('MSI1', "MSI2"), cols = c('yellow', 'red'), pt.size = 1.5, cells = 1639:1718)
FeaturePlot(gse, features = c('MSI1', "MSI2"), cols = c('yellow', 'red'), pt.size = 1.5, cells = 1719:1980)
FeaturePlot(gse, features = c('MSI1', "MSI2"), cols = c('yellow', 'red'), pt.size = 1.5, cells = 1981:2047)
FeaturePlot(gse, features = c('MSI1', "MSI2"), cols = c('yellow', 'red'), pt.size = 1.5, cells = 2048:2282)
FeaturePlot(gse, features = c('MSI1', "MSI2"), cols = c('yellow', 'red'), pt.size = 1.5, cells = 2283:2359)
FeaturePlot(gse, features = c('MSI1', "MSI2"), cols = c('yellow', 'red'), pt.size = 1.5, cells = 2360:2632)
FeaturePlot(gse, features = c('MSI1', "MSI2"), cols = c('yellow', 'red'), pt.size = 1.5, cells = 2633:2865)
FeaturePlot(gse, features = c('MSI1', "MSI2"), cols = c('yellow', 'red'), pt.size = 1.5, cells = 2866:2926)
FeaturePlot(gse, features = c('MSI1', "MSI2"), cols = c('yellow', 'red'), pt.size = 1.5, cells = 2927:3487)
FeaturePlot(gse, features = c('MSI1', "MSI2"), cols = c('yellow', 'red'), pt.size = 1.5, cells = 3488:3671)
FeaturePlot(gse, features = c('MSI1', "MSI2"), cols = c('yellow', 'red'), pt.size = 1.5, cells = 3672:4113)



# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)

# Install dittoSeq
BiocManager::install("dittoSeq", force = T)

library(dittoseq)



### Week P07 analysis: 

png(file = "UMAP_SOX2_P07.png", height = 6, width = 9, res = 1200, units = 'in')
DimPlot(object = gse, reduction = 'umap', cells = 1:170, 
        group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 07')& 
  FeaturePlot(object = gse, features = 'SOX2', cols = c('grey', 'blue'), 
              pt.size = 2, cells = 1:170, order = T, coord.fixed = T)&
  FeaturePlot(object = gse, features = 'MSI1', cols = c('grey', 'blue'), 
            pt.size = 2, cells = 1:170, order = T, coord.fixed = T)&
  FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
              pt.size = 2, cells = 1:170, order = T, coord.fixed = T)
dev.off()

png(file = "UMAP_MSI2_P07.png", height = 4.5, width = 9, res = 1200, units = 'in')
DimPlot(object = gse, reduction = 'umap', cells = 1:170, 
        group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 07')&
  FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
              pt.size = 2, cells = 1:170, order = T)
dev.off()

# Week P25 analysis: 

png(file = "UMAP_MSI1_P25.png", height = 4.5, width = 9, res = 1200, units = 'in')
DimPlot(object = gse, reduction = 'umap', cells = 3672:4113, 
        group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 25')&
  FeaturePlot(object = gse, features = 'MSI1', cols = c('grey', 'blue'), 
              pt.size = 2, cells = 3672:4113, order = T, coord.fixed = T)
dev.off()

png(file = "UMAP_MSI2_P25.png", height = 4.5, width = 9, res = 1200, units = 'in')
DimPlot(object = gse, reduction = 'umap', cells = 3672:4113, 
        group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 25')&
  FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
              pt.size = 2, cells = 3672:4113, order = T, coord.fixed = T)
dev.off()

#violin plot generation for a specific number of cells: 

# selected_cells <- names(gse$orig.ident[1:170])
# fd <- FetchData(gse,
#                   vars = 'MSI1',
#                   cells = selected_cells ,
#                   slot = "data")

#ss <- subset(x = gse[row # as features, col # as cells])

ss7 <- subset(x = gse[, 1:170])# to subset first 170 cells (weeks 7 cells only)
# as a Seurat object, which is required for vlnPlot of only 170 cells. class(ss)
ss7

png(file = "vlnPlot_MSI1_MSI2_P07_cell_type.png", height = 6, width = 8, res = 1200, units = 'in')
violin1 <- VlnPlot(object = ss7, features = 'MSI1', group.by = 'cell_type', 
        ncol = 1, pt.size =1, combine = T)& NoLegend()& labs(title = "Week 7",
                                                            subtitle = "MSI1")
violin2 <-VlnPlot(object = ss7, features = 'MSI2', group.by = 'cell_type', 
        ncol = 1, pt.size =1)& NoLegend()&labs(title="Week 7", subtitle="MSI2")

violin1&violin2

dev.off()
#https://github.com/satijalab/seurat/issues/1080



ss25 <- subset(x = gse[, 3672:4113])# to subset first 170 cells (weeks 7 cells only)
# as a Seurat object, which is required for vlnPlot of only 170 cells. class(ss)
ss25

png(file = "vlnPlot_MSI1_MSI2_P25_cell_type.png", height = 6, width = 8, res = 1200, units = 'in')
violin3 <- VlnPlot(object = ss25, features = 'MSI1', group.by = 'cell_type', 
        ncol = 1, pt.size =1, combine = T)&NoLegend()& 
  labs(title = "Week 25", subtitle = "MSI1")
violin4 <- VlnPlot(object = ss25, features = 'MSI2', group.by = 'cell_type', 
        ncol = 1, pt.size =1)&NoLegend()&
  labs(title = "Week 25", subtitle = "MSI2")

violin3&violin4
dev.off()

VlnPlot(object = ss25, features = 'MSI1') +
  stat_summary(fun = median, geom='point', size = 1, colour = "blue")



VlnPlot(object = ss25, features = 'MSI1', group.by = 'cell_type',
ncol = 1, pt.size =1, combine = T)&NoLegend()&
stat_summary(fun = mean, geom='point', size = 1, colour = "blue")

VlnPlot(object = ss25, features = 'MSI1', group.by = 'cell_type', 
        ncol = 1, pt.size =1, combine = T, )&NoLegend()&
  stat_summary(fun = rowMeans(msi1), geom='point', size = 2, colour = "blue")

mean.stat <- function(x){
  out <- rowMeans(msi1)
  names(out) <- c("ymean")
  return(out) 
}

VlnPlot(object = ss25, features = 'MSI1') +
  stat_summary(fun.y = mean.stat(), geom='point', size = 3, colour = "blue") 
# https://www.biostars.org/p/458261/
#https://github.com/satijalab/seurat/issues/2475 


###Four panel figures: 

#https://statisticsglobe.com/change-font-size-of-ggplot2-plot-in-r-axis-text-main-title-legend

#week 7: 
png(file = "UMAP_SOX2_P7.png", height = 6, width = 9, res = 1200, units = 'in')
plot1 <- DimPlot(object = gse, reduction = 'umap', cells = 1:170, 
                 group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 7')
plot2 <- FeaturePlot(object = gse, features = 'SOX2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 1:170, order = T, coord.fixed = T, combine = T)
plot3 <- FeaturePlot(object = gse, features = 'MSI1', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 1:170, order = T, coord.fixed = T, combine = T)
plot4 <- FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 1:170, order = T, coord.fixed = T, combine = T)

CombinePlots(plots = list(plot1, plot2, plot3, plot4), ncol = 2)

dev.off()

#Week 8_1
png(file = "UMAP_SOX2_P8_1.png", height = 6, width = 9, res = 1200, units = 'in')
plot1 <- DimPlot(object = gse, reduction = 'umap', cells = 171:347, 
                 group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 8_1')
plot2 <- FeaturePlot(object = gse, features = 'SOX2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 171:347, order = T, coord.fixed = T, combine = T)
plot3 <- FeaturePlot(object = gse, features = 'MSI1', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 171:347, order = T, coord.fixed = T, combine = T)
plot4 <- FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 171:347, order = T, coord.fixed = T, combine = T)

CombinePlots(plots = list(plot1, plot2, plot3, plot4), ncol = 2)

dev.off()

#week 8_2
png(file = "UMAP_SOX2_P8_2.png", height = 6, width = 9, res = 1200, units = 'in')
plot1 <- DimPlot(object = gse, reduction = 'umap', cells = 348:561, 
                 group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 8_2')
plot2 <- FeaturePlot(object = gse, features = 'SOX2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 348:561, order = T, coord.fixed = T, combine = T)
plot3 <- FeaturePlot(object = gse, features = 'MSI1', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 348:561, order = T, coord.fixed = T, combine = T)
plot4 <- FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 348:561, order = T, coord.fixed = T, combine = T)

CombinePlots(plots = list(plot1, plot2, plot3, plot4), ncol = 2)

dev.off()

#week 8_3:
png(file = "UMAP_SOX2_P8_3.png", height = 6, width = 9, res = 1200, units = 'in')
plot1 <- DimPlot(object = gse, reduction = 'umap', cells = 562:695, 
                 group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 8_3')
plot2 <- FeaturePlot(object = gse, features = 'SOX2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 562:695, order = T, coord.fixed = T, combine = T)
plot3 <- FeaturePlot(object = gse, features = 'MSI1', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 562:695, order = T, coord.fixed = T, combine = T)
plot4 <- FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 562:695, order = T, coord.fixed = T, combine = T)

CombinePlots(plots = list(plot1, plot2, plot3, plot4), ncol = 2)

dev.off()

#week 9: 
png(file = "UMAP_SOX2_P9.png", height = 6, width = 9, res = 1200, units = 'in')
plot1 <- DimPlot(object = gse, reduction = 'umap', cells = 696:866, 
                 group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 9')
plot2 <- FeaturePlot(object = gse, features = 'SOX2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 696:866, order = T, coord.fixed = T, combine = T)
plot3 <- FeaturePlot(object = gse, features = 'MSI1', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 696:866, order = T, coord.fixed = T, combine = T)
plot4 <- FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 696:866, order = T, coord.fixed = T, combine = T)

CombinePlots(plots = list(plot1, plot2, plot3, plot4), ncol = 2)

dev.off()

#week 10_1
png(file = "UMAP_SOX2_P10_1.png", height = 6, width = 9, res = 1200, units = 'in')
plot1 <- DimPlot(object = gse, reduction = 'umap', cells = 867:995, 
                 group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 10_1')
plot2 <- FeaturePlot(object = gse, features = 'SOX2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 867:995, order = T, coord.fixed = T, combine = T)
plot3 <- FeaturePlot(object = gse, features = 'MSI1', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 867:995, order = T, coord.fixed = T, combine = T)
plot4 <- FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 867:995, order = T, coord.fixed = T, combine = T)

CombinePlots(plots = list(plot1, plot2, plot3, plot4), ncol = 2)

dev.off()

#week 10_2
png(file = "UMAP_SOX2_P10_2.png", height = 6, width = 9, res = 1200, units = 'in')
plot1 <- DimPlot(object = gse, reduction = 'umap', cells = 996:1231, 
                 group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 10_2')
plot2 <- FeaturePlot(object = gse, features = 'SOX2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 996:1231, order = T, coord.fixed = T, combine = T)
plot3 <- FeaturePlot(object = gse, features = 'MSI1', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 996:1231, order = T, coord.fixed = T, combine = T)
plot4 <- FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 996:1231, order = T, coord.fixed = T, combine = T)

CombinePlots(plots = list(plot1, plot2, plot3, plot4), ncol = 2)

dev.off()

#week 13
png(file = "UMAP_SOX2_P13.png", height = 6, width = 9, res = 1200, units = 'in')
plot1 <- DimPlot(object = gse, reduction = 'umap', cells = 1232:1384, 
                 group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 13')
plot2 <- FeaturePlot(object = gse, features = 'SOX2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 1232:1384, order = T, coord.fixed = T, combine = T)
plot3 <- FeaturePlot(object = gse, features = 'MSI1', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 1232:1384, order = T, coord.fixed = T, combine = T)
plot4 <- FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 1232:1384, order = T, coord.fixed = T, combine = T)

CombinePlots(plots = list(plot1, plot2, plot3, plot4), ncol = 2)

dev.off()

#week 14
png(file = "UMAP_SOX2_P14.png", height = 6, width = 9, res = 1200, units = 'in')
plot1 <- DimPlot(object = gse, reduction = 'umap', cells = 1385:1563, 
                 group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 14')
plot2 <- FeaturePlot(object = gse, features = 'SOX2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 1385:1563, order = T, coord.fixed = T, combine = T)
plot3 <- FeaturePlot(object = gse, features = 'MSI1', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 1385:1563, order = T, coord.fixed = T, combine = T)
plot4 <- FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 1385:1563, order = T, coord.fixed = T, combine = T)

CombinePlots(plots = list(plot1, plot2, plot3, plot4), ncol = 2)

dev.off()

#week 15
png(file = "UMAP_SOX2_P15.png", height = 6, width = 9, res = 1200, units = 'in')
plot1 <- DimPlot(object = gse, reduction = 'umap', cells = 1564:1638, 
                 group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 15')
plot2 <- FeaturePlot(object = gse, features = 'SOX2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 1564:1638, order = T, coord.fixed = T, combine = T)
plot3 <- FeaturePlot(object = gse, features = 'MSI1', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 1564:1638, order = T, coord.fixed = T, combine = T)
plot4 <- FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 1564:1638, order = T, coord.fixed = T, combine = T)

CombinePlots(plots = list(plot1, plot2, plot3, plot4), ncol = 2)

dev.off()

#week 16
png(file = "UMAP_SOX2_P16.png", height = 6, width = 9, res = 1200, units = 'in')
plot1 <- DimPlot(object = gse, reduction = 'umap', cells = 1639:1718, 
                 group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 16')
plot2 <- FeaturePlot(object = gse, features = 'SOX2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 1639:1718, order = T, coord.fixed = T, combine = T)
plot3 <- FeaturePlot(object = gse, features = 'MSI1', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 1639:1718, order = T, coord.fixed = T, combine = T)
plot4 <- FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 1639:1718, order = T, coord.fixed = T, combine = T)

CombinePlots(plots = list(plot1, plot2, plot3, plot4), ncol = 2)

dev.off()

#week 17: 
png(file = "UMAP_SOX2_P17.png", height = 6, width = 9, res = 1200, units = 'in')
plot1 <- DimPlot(object = gse, reduction = 'umap', cells = 1719:1980, 
                 group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 17')
plot2 <- FeaturePlot(object = gse, features = 'SOX2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 1719:1980, order = T, coord.fixed = T, combine = T)
plot3 <- FeaturePlot(object = gse, features = 'MSI1', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 1719:1980, order = T, coord.fixed = T, combine = T)
plot4 <- FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 1719:1980, order = T, coord.fixed = T, combine = T)

CombinePlots(plots = list(plot1, plot2, plot3, plot4), ncol = 2)

dev.off()

#week 19_1
png(file = "UMAP_SOX2_P19_1.png", height = 6, width = 9, res = 1200, units = 'in')
plot1 <- DimPlot(object = gse, reduction = 'umap', cells = 1981:2047, 
                 group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 19_1')
plot2 <- FeaturePlot(object = gse, features = 'SOX2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 1981:2047, order = T, coord.fixed = T, combine = T)
plot3 <- FeaturePlot(object = gse, features = 'MSI1', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 1981:2047, order = T, coord.fixed = T, combine = T)
plot4 <- FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 1981:2047, order = T, coord.fixed = T, combine = T)

CombinePlots(plots = list(plot1, plot2, plot3, plot4), ncol = 2)

dev.off()

#week 19_2
png(file = "UMAP_SOX2_P19_2.png", height = 6, width = 9, res = 1200, units = 'in')
plot1 <- DimPlot(object = gse, reduction = 'umap', cells = 2048:2282, 
                 group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 19_2')
plot2 <- FeaturePlot(object = gse, features = 'SOX2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 2048:2282, order = T, coord.fixed = T, combine = T)
plot3 <- FeaturePlot(object = gse, features = 'MSI1', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 2048:2282, order = T, coord.fixed = T, combine = T)
plot4 <- FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 2048:2282, order = T, coord.fixed = T, combine = T)

CombinePlots(plots = list(plot1, plot2, plot3, plot4), ncol = 2)

dev.off()

#week 21
png(file = "UMAP_SOX2_P21.png", height = 6, width = 9, res = 1200, units = 'in')
plot1 <- DimPlot(object = gse, reduction = 'umap', cells = 2283:2359, 
                 group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 21')
plot2 <- FeaturePlot(object = gse, features = 'SOX2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 2283:2359, order = T, coord.fixed = T, combine = T)
plot3 <- FeaturePlot(object = gse, features = 'MSI1', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 2283:2359, order = T, coord.fixed = T, combine = T)
plot4 <- FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 2283:2359, order = T, coord.fixed = T, combine = T)

CombinePlots(plots = list(plot1, plot2, plot3, plot4), ncol = 2)

dev.off()

#week 22_1
png(file = "UMAP_SOX2_P22_1.png", height = 6, width = 9, res = 1200, units = 'in')
plot1 <- DimPlot(object = gse, reduction = 'umap', cells = 2360:2632, 
                 group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 22_1')
plot2 <- FeaturePlot(object = gse, features = 'SOX2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 2360:2632, order = T, coord.fixed = T, combine = T)
plot3 <- FeaturePlot(object = gse, features = 'MSI1', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 2360:2632, order = T, coord.fixed = T, combine = T)
plot4 <- FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 2360:2632, order = T, coord.fixed = T, combine = T)

CombinePlots(plots = list(plot1, plot2, plot3, plot4), ncol = 2)

dev.off()

#week 22_2
png(file = "UMAP_SOX2_P22_2.png", height = 6, width = 9, res = 1200, units = 'in')
plot1 <- DimPlot(object = gse, reduction = 'umap', cells = 2633:2865, 
                 group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 22_2')
plot2 <- FeaturePlot(object = gse, features = 'SOX2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 2633:2865, order = T, coord.fixed = T, combine = T)
plot3 <- FeaturePlot(object = gse, features = 'MSI1', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 2633:2865, order = T, coord.fixed = T, combine = T)
plot4 <- FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 2633:2865, order = T, coord.fixed = T, combine = T)

CombinePlots(plots = list(plot1, plot2, plot3, plot4), ncol = 2)

dev.off()

#week 23_1
png(file = "UMAP_SOX2_P23_1.png", height = 6, width = 9, res = 1200, units = 'in')
plot1 <- DimPlot(object = gse, reduction = 'umap', cells = 2866:2926, 
                 group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 23_1')
plot2 <- FeaturePlot(object = gse, features = 'SOX2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 2866:2926, order = T, coord.fixed = T, combine = T)
plot3 <- FeaturePlot(object = gse, features = 'MSI1', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 2866:2926, order = T, coord.fixed = T, combine = T)
plot4 <- FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 2866:2926, order = T, coord.fixed = T, combine = T)

CombinePlots(plots = list(plot1, plot2, plot3, plot4), ncol = 2)

dev.off()

#week 23_2
png(file = "UMAP_SOX2_P23_2.png", height = 6, width = 9, res = 1200, units = 'in')
plot1 <- DimPlot(object = gse, reduction = 'umap', cells = 2927:3487, 
                 group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 23_2')
plot2 <- FeaturePlot(object = gse, features = 'SOX2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 2927:3487, order = T, coord.fixed = T, combine = T)
plot3 <- FeaturePlot(object = gse, features = 'MSI1', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 2927:3487, order = T, coord.fixed = T, combine = T)
plot4 <- FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 2927:3487, order = T, coord.fixed = T, combine = T)

CombinePlots(plots = list(plot1, plot2, plot3, plot4), ncol = 2)

dev.off()

#week 24: 
png(file = "UMAP_SOX2_P24.png", height = 6, width = 9, res = 1200, units = 'in')
plot1 <- DimPlot(object = gse, reduction = 'umap', cells = 3488:3671, 
                 group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 24')
plot2 <- FeaturePlot(object = gse, features = 'SOX2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 3488:3671, order = T, coord.fixed = T, combine = T)
plot3 <- FeaturePlot(object = gse, features = 'MSI1', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 3488:3671, order = T, coord.fixed = T, combine = T)
plot4 <- FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 3488:3671, order = T, coord.fixed = T, combine = T)

CombinePlots(plots = list(plot1, plot2, plot3, plot4), ncol = 2)

dev.off()

#week 25:
png(file = "UMAP_SOX2_P25.png", height = 6, width = 9, res = 1200, units = 'in')
plot1 <- DimPlot(object = gse, reduction = 'umap', cells = 3672:4113, 
                 group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 25')
plot2 <- FeaturePlot(object = gse, features = 'SOX2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 3672:4113, order = T, coord.fixed = T, combine = T)
plot3 <- FeaturePlot(object = gse, features = 'MSI1', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 3672:4113, order = T, coord.fixed = T, combine = T)
plot4 <- FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 3672:4113, order = T, coord.fixed = T, combine = T)

CombinePlots(plots = list(plot1, plot2, plot3, plot4), ncol = 2)

dev.off()

#male only: 
male <- c(171:347, 562:695, 696:866, 1232:1384, 1719:1980, 2283:2359, 2360:2632, 
          2866:2926, 2927:3487, 3672:4113)
female <- c(1:170, 348:561, 867:995, 996:1231, 1385:1563, 1564:1638, 1639:1718, 
            1981:2047, 2048:2282, 2633:2865, 3488:3671)

DimPlot(object = gse, reduction = 'umap', cells = female, 
        group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Female')

DimPlot(object = gse, reduction = 'umap', pt.size = 2, cells.highlight = list(male, female),
        cols.highlight = c('skyblue', 'magenta'), order = F, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Female')

view(gse@meta.data)

#https://github.com/satijalab/seurat/issues/2371
#g1_untreat <- WhichCells(integrated, idents = c("group1_untreated")
#g1_treat <- WhichCells(integrated, idents = c( "group1_treated")
#DimPlot(integrated, label=T, group.by="Treat", cells.highlight= list(g1_treat, g1_untreat), cols.highlight = c("darkblue", "darkred"), cols= "grey")


#Violin plots for cells subsets: week 7 to week 25:

ss7 <- subset(x = gse[, 1:170])
ss8_1 <- subset(x = gse[, 171:347])
ss8_2 <- subset(x = gse[, 348:561])
ss8_3 <- subset(x = gse[, 562:695])
ss9 <- subset(x = gse[, 696:866])
ss10_1 <- subset(x = gse[, 867:995])
ss10_2 <- subset(x = gse[, 996:1231])
ss13 <- subset(x = gse[, 1232:1384])
ss14 <- subset(x = gse[, 1385:1563])
ss15 <- subset(x = gse[, 1564:1638])
ss16 <- subset(x = gse[, 1639:1718])
ss17 <- subset(x = gse[, 1719:1980])
ss19_1 <- subset(x = gse[, 1981:2047])
ss19_2 <- subset(x = gse[, 2048:2282])
ss21 <- subset(x = gse[, 2283:2359])
ss22_1 <- subset(x = gse[, 2360:2632])
ss22_2 <- subset(x = gse[, 2633:2865])
ss23_1<- subset(x = gse[, 2866:2926])
ss23_2 <- subset(x = gse[, 2927:3487])
ss24 <- subset(x = gse[, 3488:3671])
ss25 <- subset(x = gse[, 3672:4113])

#https://stackoverflow.com/questions/8399100/r-plot-size-and-resolution

#Week 7: 
png(file = "VlnPlot_SOX2_P7_cell_type.png", height = 9, width = 11, res = 1200, units = 'in')
violin1 <- VlnPlot(object = ss7, features = 'SOX2', group.by = 'cell_type', 
                   ncol = 1, pt.size =1, combine = T, same.y.lims = T)& 
  NoLegend()& labs(title = "Week 7", subtitle = "SOX2")

violin2 <-VlnPlot(object = ss7, features = 'MSI1', group.by = 'cell_type', 
                  ncol = 1, pt.size =1, same.y.lims = T)& 
  NoLegend()&labs(title=" ", subtitle="MSI1")

violin3 <-VlnPlot(object = ss7, features = 'MSI2', group.by = 'cell_type', 
                  ncol = 1, pt.size =1, same.y.lims = T)& 
  NoLegend()&labs(title=" ", subtitle="MSI2")

violin1+violin2+violin3

dev.off()

#Week 8_1: 
png(file = "VlnPlot_SOX2_P8_1_cell_type.png", height = 9, width = 11, res = 1200, units = 'in')
violin1 <- VlnPlot(object = ss8_1, features = 'SOX2', group.by = 'cell_type', 
                   ncol = 1, pt.size =1, combine = T)& NoLegend()& labs(title = "Week 8_1",
                                                                        subtitle = "SOX2")
violin2 <-VlnPlot(object = ss8_1, features = 'MSI1', group.by = 'cell_type', 
                  ncol = 1, pt.size =1)& NoLegend()&labs(title=" ", subtitle="MSI1")

violin3 <-VlnPlot(object = ss8_1, features = 'MSI2', group.by = 'cell_type', 
                  ncol = 1, pt.size =1)& NoLegend()&labs(title=" ", subtitle="MSI2")

violin1+violin2+violin3

dev.off()

#Week 8_2: 
png(file = "VlnPlot_SOX2_P8_2_cell_type.png", height = 9, width = 11, res = 1200, units = 'in')
violin1 <- VlnPlot(object = ss8_2, features = 'SOX2', group.by = 'cell_type', 
                   ncol = 1, pt.size =1, combine = T)& NoLegend()& labs(title = "Week 8_2",
                                                                        subtitle = "SOX2")
violin2 <-VlnPlot(object = ss8_2, features = 'MSI1', group.by = 'cell_type', 
                  ncol = 1, pt.size =1)& NoLegend()&labs(title=" ", subtitle="MSI1")

violin3 <-VlnPlot(object = ss8_2, features = 'MSI2', group.by = 'cell_type', 
                  ncol = 1, pt.size =1)& NoLegend()&labs(title=" ", subtitle="MSI2")

violin1+violin2+violin3

dev.off()

#Week 8_3: 
#https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/
png(file = "VlnPlot_SOX2_P8_3_cell_type.png", height = 9, width = 11, res = 1200, units = 'in')
violin1 <- VlnPlot(object = ss8_3, features = c('SOX2', 'MSI1', 'MSI2'), group.by = 'cell_type',
                   pt.size =1, same.y.lims = T,stack = T, flip = T)& 
  NoLegend()& labs(title = "Week 8_3")

violin2 <-VlnPlot(object = ss8_3, features = 'MSI1', group.by = 'cell_type', 
                  ncol = 1, pt.size =1, same.y.lims = T)& 
  NoLegend()&labs(title=" ", subtitle="MSI1")

violin3 <-VlnPlot(object = ss8_3, features = 'MSI2', group.by = 'cell_type', 
                  ncol = 1, pt.size =1, same.y.lims = T)& 
  NoLegend()&labs(title=" ", subtitle="MSI2")

violin1+violin2+violin3

dev.off()


VlnPlot(object = ss8_3, features = c('SOX2', 'MSI1', 'MSI2'), group.by = 'cell_type',
        pt.size =1, same.y.lims = T,stack = T, flip = T)+labs(title = "Week 8_3")+
  theme(axis.text = element_text(size = 5))



###Five panel figures: 
#week 7: 
png(file = "UMAP_SOX2_PROP1_P7.png", height = 8, width = 14, res = 1200, units = 'in')
plot1 <- DimPlot(object = gse, reduction = 'umap', cells = 1:170, 
                 group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 7')
plot2 <- FeaturePlot(object = gse, features = 'SOX2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 1:170, order = T, coord.fixed = T, combine = T)
plot3 <-   FeaturePlot(object = gse, features = 'PROP1', cols = c('grey', 'blue'), 
                       pt.size = 2, cells = 1:170, order = T, coord.fixed = T, combine = T)
plot4 <- FeaturePlot(object = gse, features = 'MSI1', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 1:170, order = T, coord.fixed = T, combine = T)
plot5 <- FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 1:170, order = T, coord.fixed = T, combine = T)

CombinePlots(plots = list(plot1, plot2, plot3, plot4, plot5), ncol = 3)

dev.off()

#week 25: 
png(file = "UMAP_SOX2_PROP1_P25.png", height = 8, width = 14, res = 1200, units = 'in')
plot1 <- DimPlot(object = gse, reduction = 'umap', cells = 3672:4113, 
                 group.by = 'cell_type', pt.size = 2, order = T, combine = T)&
  coord_fixed(ratio = 1)& labs(title = 'Week 25')
plot2 <- FeaturePlot(object = gse, features = 'SOX2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 3672:4113, order = T, coord.fixed = T, combine = T)
plot3 <-   FeaturePlot(object = gse, features = 'PROP1', cols = c('grey', 'blue'), 
                       pt.size = 2, cells = 3672:4113, order = T, coord.fixed = T, combine = T)
plot4 <- FeaturePlot(object = gse, features = 'MSI1', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 3672:4113, order = T, coord.fixed = T, combine = T)
plot5 <- FeaturePlot(object = gse, features = 'MSI2', cols = c('grey', 'blue'), 
                     pt.size = 2, cells = 3672:4113, order = T, coord.fixed = T, combine = T)

CombinePlots(plots = list(plot1, plot2, plot3, plot4, plot5), ncol = 3)

dev.off()


#violin plots: 
ss7 <- subset(x = gse[, 1:170])

#Week 7: 
png(file = "VlnPlot_SOX2_PROP1_P7_cell_type.png", height = 12, width = 12, res = 1200, units = 'in')
violin1 <- VlnPlot(object = ss7, features = 'SOX2', group.by = 'cell_type', 
                   ncol = 1, pt.size =1.5, combine = T, same.y.lims = T)& 
  NoLegend()& labs(title = "Week 7", subtitle = "SOX2")

violin2 <- VlnPlot(object = ss7, features = 'PROP1', group.by = 'cell_type', 
                   ncol = 1, pt.size =1.5, combine = T, same.y.lims = T)& 
  NoLegend()& labs(title = " ", subtitle = "PROP1")

violin3 <-VlnPlot(object = ss7, features = 'MSI1', group.by = 'cell_type', 
                  ncol = 1, pt.size =1.5, same.y.lims = T)& 
  NoLegend()&labs(title=" ", subtitle="MSI1")

violin4 <-VlnPlot(object = ss7, features = 'MSI2', group.by = 'cell_type', 
                  ncol = 1, pt.size =1.5, same.y.lims = T)& 
  NoLegend()&labs(title=" ", subtitle="MSI2")

violin1+violin2+violin3+violin4

dev.off()

#week 25
ss25 <- subset(x = gse[, 3672:4113])

png(file = "VlnPlot_SOX2_PROP1_P25_cell_type.png", height = 12, width = 16, res = 1200, units = 'in')
violin1 <- VlnPlot(object = ss25, features = 'SOX2', group.by = 'cell_type', 
                   ncol = 1, pt.size =1.5, combine = T, same.y.lims = T)& 
  NoLegend()& labs(title = "Week 25", subtitle = "SOX2")

violin2 <- VlnPlot(object = ss25, features = 'PROP1', group.by = 'cell_type', 
                   ncol = 1, pt.size =1.5, combine = T, same.y.lims = T)& 
  NoLegend()& labs(title = " ", subtitle = "PROP1")

violin3 <-VlnPlot(object = ss25, features = 'MSI1', group.by = 'cell_type', 
                  ncol = 1, pt.size =1.5, same.y.lims = T)& 
  NoLegend()&labs(title=" ", subtitle="MSI1")

violin4 <-VlnPlot(object = ss25, features = 'MSI2', group.by = 'cell_type', 
                  ncol = 1, pt.size =1.5, same.y.lims = T)& 
  NoLegend()&labs(title=" ", subtitle="MSI2")

violin1+violin2+violin3+violin4

dev.off()


#DotPlot: 
#week 25: 
ss25 <- subset(x = gse[, 3672:4113])

png(file = "DotPlot_cell_types_P25.png", height = 5, width = 12, res = 1200, units = 'in')
DotPlot(ss25, features = c("SOX2","PROP1", "POU1F1", 'MSI1', 'MSI2', 'POMC', 
                           'GATA2', 'NR5A1', 'TSHB', 'PRL', 'FSHB', 'GNRHR'), 
        group.by = "cell_type")+labs(title = "Week 25")+theme_bw()
dev.off()

png(file = "DotPlot_cell_types_P25_5percentFiltered.png", height = 5, width = 12, res = 1200, units = 'in')
DotPlot(ss25, features = c("SOX2","PROP1", "POU1F1", 'MSI1', 'MSI2', 'POMC', 
                           'GATA2', 'NR5A1', 'TSHB', 'PRL', 'FSHB', 'GNRHR'), 
        group.by = "cell_type", dot.min = 0.05)+labs(title = "Week 25")+theme_bw()
dev.off()

png(file = "DotPlot_cell_types_P25_5percentFiltered_NoPRL.png", height = 5, width = 10, res = 1200, units = 'in')
DotPlot(ss25, features = c("SOX2","PROP1", "POU1F1", 'MSI1', 'MSI2', 'POMC', 
                           'GATA2', 'NR5A1', 'TSHB', 'FSHB', 'GNRHR'), 
        group.by = "cell_type", dot.min = 0.05)+labs(title = "Week 25")+theme_bw()
dev.off()


###MSI1 and MSI2 subsets: 
mean(ss7@assays$RNA['MSI1']) #0.8523797
mean(ss7@assays$RNA['MSI2']) #0.4097002
mean(ss25@assays$RNA['MSI1']) #0.2143217
mean(ss25@assays$RNA['MSI2']) #0.264627

length(which(ss7@assays$RNA['MSI1']>0)) #105
length(which(ss7@assays$RNA['MSI2']>0)) #49
length(which(ss25@assays$RNA['MSI1']>0)) #69
length(which(ss25@assays$RNA['MSI2']>0)) #73

#subsetting MSI1 and MSI2 expressing cells only in week 7 & week 25: 
ss7msi1 <- subset(x = ss7[, which(ss7@assays$RNA['MSI1']>0)])
ss7msi2 <- subset(x = ss7[, which(ss7@assays$RNA['MSI2']>0)])
ss25msi1 <- subset(x = ss25[, which(ss25@assays$RNA['MSI1']>0)])
ss25msi2 <- subset(x = ss25[, which(ss25@assays$RNA['MSI2']>0)])

mean(ss7msi1@assays$RNA['MSI1']) #1.380043
mean(ss7msi2@assays$RNA['MSI2']) #1.421409
mean(ss25msi1@assays$RNA['MSI1']) #1.372902
mean(ss25msi2@assays$RNA['MSI2']) #1.602262

df1 <- t(as.data.frame(ss7msi1@assays$RNA['MSI1']))
df2 <- t(as.data.frame(ss7msi2@assays$RNA['MSI2']))
df3 <- t(as.data.frame(ss25msi1@assays$RNA['MSI1']))
df4 <- t(as.data.frame(ss25msi2@assays$RNA['MSI2']))

dfAll <- rbind(df1, df2, df3, df4)
colnames(dfAll) <- NULL

msi1 <- t(as.data.frame(data@assays$RNA['MSI1']))
msi2 <- t(as.data.frame(data@assays$RNA['MSI2']))
msi <- cbind(msi1, msi2)
msi170 <- msi[1:170, ] #for week 7
msi3672 <- msi[3672:4113, ] #for week 25

length(which(msi170[, 'MSI1']>msi170[, 'MSI2'])) #85cells
length(which(msi170[, 'MSI1']<msi170[, 'MSI2'])) #27 cells
length(which(msi170[, 'MSI1']==msi170[, 'MSI2'])) #58 cells

length(which(msi3672[, 'MSI1']>msi3672[, 'MSI2'])) #53cells
length(which(msi3672[, 'MSI1']<msi3672[, 'MSI2'])) #70 cells
length(which(msi3672[, 'MSI1']==msi3672[, 'MSI2'])) #319 cells

AverageExpression(object = ss7, features = c('MSI1', 'MSI2'), 
                  group.by = 'cell_type', slot = 'data') #calculates avg. expression
#on the logNormalized data
AverageExpression(object = ss7, features = c('MSI1', 'MSI2'), 
                  group.by = 'cell_type', slot = 'counts')##calculates avg. expression
#on the non-logNormalized data

AverageExpression(object = ss25, features = c('MSI1', 'MSI2'), 
                  group.by = 'cell_type', slot = 'data')


###CellChat: https://rdrr.io/github/sqjin/CellChat/f/tutorial/CellChat-vignette.Rmd
install.packages('remotes')
BiocManager::install("ComplexHeatmap")
library(remotes)
library(ComplexHeatmap)
remotes::install_github('sqjin/CellChat')
library(CellChat)
options(StringsAsFactors=T)

# Stash cell identity classes
gse[["Idents4CellChat"]] <- Idents(object = gse) #CreateCellChat function doesn't accept
#cell labels as 0. So, the plan is I'm copying seurat_clusters from the metadata
#and creating a new column names 'Idents4CellChat' where I will rename label 0 to
# 21 and use this column for group.by argument.


length(which(gse$Idents4CellChat==0))#to figure out how many cells with factor 0 in that column
gse$Idents4CellChat <- NULL #to delete newly created metadata column. 

# Rename identity classes
gse2 <- gse
gse2 <- RenameIdents(object = gse2, `0` = "21")

#Create a CellChat object: 
CellChat <- createCellChat(object = gse, group.by = 'cell_type', assay = 'RNA')

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
CellChat@DB <- CellChatDB #the whole database

# subset the expression data of signaling genes for saving computation cost
CellChat <- subsetData(CellChat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
CellChat <- identifyOverExpressedGenes(CellChat)
CellChat <- identifyOverExpressedInteractions(CellChat)
# project gene expression data onto PPI network (optional)
CellChat <- projectData(CellChat, PPI.human) #useful when working with shallow-depth 
#sequencing data. 

CellChat <- computeCommunProb(CellChat, raw.use = T, population.size = T) #raw.use=T when 
#the sequencing depth is higher. 
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
CellChat <- filterCommunication(CellChat, min.cells = 10)

# Extract the inferred cellular communication network as a data frame
# 
# We provide a function subsetCommunication to easily access the inferred cell-cell 
#communications of interest. For example,
# 
# df.net <- subsetCommunication(cellchat) returns a data frame consisting of all 
#the inferred cell-cell communications at the level of ligands/receptors. 
#Set slot.name = "netP" to access the the inferred communications at the level of 
#signaling pathways
# 
# df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) 
#gives the inferred cell-cell communications sending from cell groups 1 and 2 to 
#cell groups 4 and 5.
# 
# df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) gives 
#the inferred cell-cell communications mediated by signaling WNT and TGFb.

CellChat <- computeCommunProbPathway(CellChat)
CellChat <- aggregateNet(CellChat)

groupSize <- as.numeric(table(CellChat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(CellChat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(CellChat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- CellChat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

###Visualization: 

pathways.show <- c("CXCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(CellChat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(CellChat, signaling = pathways.show, layout = "circle")
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(CellChat, signaling = pathways.show, layout = "chord")
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(CellChat, signaling = pathways.show, color.heatmap = "Reds")

netAnalysis_contribution(CellChat, signaling = pathways.show)

#Automatically save the plots of the all inferred network for quick exploration

#In practical use, USERS can use 'for ... loop' to automatically save the all inferred network 
#for quick exploration using netVisual. netVisual supports an output in the formats of svg, png and pdf.

# Access all the signaling pathways showing significant communications
pathways.show.all <- CellChat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(CellChat@idents)
vertex.receiver = seq(1, 4)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(CellChat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(CellChat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

netVisual_bubble(CellChat, sources.use = 13, targets.use = c(5:11), remove.isolate = FALSE)

netVisual_bubble(CellChat, sources.use =13, targets.use = c(5:11), signaling = c("MK", "NOTCH", "GH"), remove.isolate = FALSE)

plotGeneExpression(CellChat, signaling = "CXCL")


# Plot order VlnPlot #454
# How can I change the x-axis ordering in VlnPlot? #1250
###to have violin plot with sorted cell types: load the .rds file first. 
# Define an order of cluster identities
levels(as.factor(gse@meta.data$cell_type))
my_levels <- c('Somatotrope', 'Lactotrope', 'Thyrotrope', 'Corticotrope', 
               'Gonadotrope', 'Pre.Gonado', 'Pro.PIT1', 'Stem', 'CC', 
               'PP', 'MC', 'EC', 'Imm', 'RBC')

# Relevel object@ident
gse@meta.data$cell_type <- factor(x = gse@meta.data$cell_type,
                                         levels = my_levels)


features <- c("MSI1","MSI2", "ORMDL1", 'ORMDL3', 'INSIG1', 'RAB2A', 
              'RAB3D', 'RAB21', 'RAB27A', 'ATF4', 'CREB3L2')
png(file = "VlnPlot_cell_types_sorted.png", height = 7, width = 10, res = 1200, units = 'in')
violin1 <- VlnPlot(object = gse, features = features, pt.size =1.5, group.by = 'cell_type',
                   combine = T, same.y.lims = T, stack = T, flip = T)& 
  NoLegend()& labs(title = "Human Fetal Pituitary")

violin1
dev.off()

###NORAD aka LINC00657 https://www.ncbi.nlm.nih.gov/gene/647979 

png(file = "VlnPlot_cell_types_NORAD.png", height = 6, width = 8, res = 1200, units = 'in')
NORAD <- VlnPlot(object = gse, features = c("MSI1", "MSI2", "LINC00657"), pt.size =1.5, group.by = 'cell_type',
                   combine = T, same.y.lims = T, stack = T, flip = T)& NoLegend()& labs(title = "Human Fetal Pituitary")

NORAD
dev.off()

FeaturePlot(object = gse, features = 'LINC00657', min.cutoff = 'q10')

png(file = "DotPlot_cell_types_NORAD.png", height = 5, width = 6, res = 1200, units = 'in')
DotPlot(gse, features = c("MSI1","MSI2", 'LINC00657'), cols = c('grey', 'red'),
        group.by = "cell_type", dot.min = 0.05)+labs(title = "Human Fetal Pituitary")+theme_bw()
dev.off()

#https://bioinformatics.stackexchange.com/questions/10811/seurat-vlnplot-presenting-expression-of-multiple-genes-in-a-single-cluster
#this thread describes vlnPlot with multiple features in the same cells types. 

###this thread describes how to find out percent of gene expressing cells 
# https://github.com/satijalab/seurat/issues/371


