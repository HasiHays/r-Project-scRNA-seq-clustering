#seurat scRNA-seq working flow
#setwd("~/Desktop/R/Seurat-workflow/script")
#data - NSCLC data form 10X Genomics

#(A)- Pre processing and clusteering
# 01. QC-------------------------------------------
# 02. Filtering -------------------------------------------
# 03. Normalize-------------------------------------------
# 04. Identify highly variable features-------------------------------------
# 5. Scaling---------------------------
# 6. Linear dimensional reduction-------
#7. clustering ---------------
#####
#Number of reads---------------------
#########


#(B)- Analyse using Seurat Object

#Take libraries
library(Seurat)
library(tidyverse)

#load data - NSCLC data form 10X Genomics
my.data <- Read10X_h5(filename = "../data/40k_NSCLC_DTC_3p_HT_nextgem_Multiplex_count_raw_feature_bc_matrix.h5")
cts <- my.data$`Gene Expression`

#make seurat 
my.data.obj <-CreateSeuratObject(counts = cts, project = "seuratProject", min.cells = 3, min.features = 200)


# 01. QC-------------------------------------------
View(my.data.obj@meta.data)
# % of MT reads
### check MT genes
grep("^MT-",rownames(my.data.obj@assays$RNA@counts),value = TRUE)

my.data.obj$present.MT <- PercentageFeatureSet(my.data.obj, pattern = "^MT-")
head(my.data.obj)

VlnPlot(my.data.obj, features = c("nFeature_RNA", "nCount_RNA", "present.MT"), ncol = 3)

####Check ribosome
grep("^RP[LS]",rownames(my.data.obj@assays$RNA@counts),value = TRUE)

my.data.obj$present.Ribosomal <- PercentageFeatureSet(my.data.obj, pattern = "^RP[LS]")
head(my.data.obj)

VlnPlot(my.data.obj, features = c("nFeature_RNA", "nCount_RNA", "present.MT", "present.Ribosomal" ), ncol = 4)

####
FeatureScatter(my.data.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')



# 02. Filtering -------------------------------------------
my.data.obj <- subset(my.data.obj, subset = nFeature_RNA >200 &  nFeature_RNA >2500 & present.MT < 5)
### if want, filter ribosome and doublets too


# 03. Normalize-------------------------------------------
my.data.obj <- NormalizeData(my.data.obj,
                             normalization.method = "LogNormalize", scale.factor = 100000)

### or -
my.data.obj <- NormalizeData(my.data.obj)


# 04. Identify highly variable features-------------------------------------
my.data.obj <- FindVariableFeatures(my.data.obj,selection.method = "vst", nfeatures = 2000)

##top 10 genes
top10 <- head(VariableFeatures(my.data.obj), 10)

plot1 <- VariableFeaturePlot(my.data.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

# 5. Scaling---------------------------
all.genes <- rownames(my.data.obj)
my.data.obj <- ScaleData(my.data.obj, features = all.genes)

str(my.data.obj)

# 6. Linear dimensional reduction-------PCA=Principle component analysis-identify hetrogenety
my.data.obj <- RunPCA(my.data.obj, features = VariableFeatures(object = my.data.obj))

## visualize PCA results -top 5
print(my.data.obj[["pca"]], dims = 1:15, nfeatures = 5)
DimHeatmap(my.data.obj, dims = 1, cells = 500, balanced = TRUE)

## determine dimentionality of the data 
ElbowPlot(my.data.obj) 


#7. clustering ---------------
my.data.obj <- FindNeighbors(my.data.obj, dims = 1:15)

#understand resolution ----- high res- more cluster
my.data.obj <- FindClusters(my.data.obj, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
view(my.data.obj@meta.data)

DimPlot(my.data.obj, group.by = "RNA_snn_res.0.5", label = TRUE)

## setting identity of the clusters
Idents(my.data.obj)
Idents(my.data.obj) <- "RNA_snn_res.0.1"
Idents(my.data.obj)

DimPlot(my.data.obj, group.by = "RNA_snn_res.0.1", label = TRUE)
##non linear dimensionality reduction
my.data.obj <-RunUMAP(my.data.obj, dims = 1:15)

DimPlot(my.data.obj, reduction = "umap")

#save file t for furter analysis
saveRDS(my.data.obj, file = "../data/my.data_processed.rds")

################################################################
#Number of reads---------------------

head(my.data.obj)
VlnPlot(my.data.obj, features = "nCount_RNA")
VlnPlot(my.data.obj, features = "nFeature_RNA")
VlnPlot(my.data.obj, features = "present.MT")

##########################################
# (B)- Analyse using Seurat Object
############################################
# 1. 
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
My.data.final <- readRDS("../data/my.data_processed.rds")

My.data.final$groups <- sample(c("group1", "group2"), size = ncol(My.data.final), replace = TRUE)
features <- c( "IGKC" ,   "SFTPC" ,  "IGHG1" ,  "IGHA1"  , "SCGB1A1", "SCGB3A1")
My.data.final

# Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
RidgePlot(My.data.final, features = features, ncol = 2)

# Violin plot - Visualize single cell expression distributions in each cluster
VlnPlot(My.data.final, features = features)

# Feature plot - visualize feature expression in low-dimensional space
FeaturePlot(My.data.final, features = features)

# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
DotPlot(My.data.final, features = features) + RotatedAxis()

# Single cell heatmap of feature expression
DoHeatmap(subset(My.data.final, downsample = 100), features = features, size = 3)

# Plot a legend to map colors to expression levels
FeaturePlot(My.data.final, features = "IGKC")

# Adjust the contrast in the plot
FeaturePlot(My.data.final, features = "IGKC", min.cutoff = 1, max.cutoff = 3)





