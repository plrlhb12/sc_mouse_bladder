# 1. the starting point is the standard unsuperivised clustering workflow
it includes slection and filtration of cells based on QC metric, data normalization and scaling, and the detection of highly variable features

# Notice the terminologies using by Seurat

# **multidmodal analysis**, simultaneous measurements of multiple data types from the same cell
# e.g., CITE-seq has both GEX and antibody barcode; 10x multiome kit has both GEX and ATAC

# **Data integration** means integrate analysis generated acrooss differnt conditon, tech, or species which contain batch effects
# it has several methods for integration, 
e.g., 
1. Find anchor
2. mapping and annotating query dataset onto a reference, 
3. indentify anchors ssuing reciprocal PCA (rPCA, faster and more conderved)
4. BPCells for integrate very large datasets (> 200,000 cells)
5. integrate scRNA-Sdeq and scATAC dataset
6. multimodal reference atlases

# 2. more new methods
1. WNN: analyze multimodal sc data with weighted nearest neighbor analysis, e.g, define cell state based on multiple modalities
2. Mixscape: for pooled sc CRISPR perturbation screen
3. SCTransform: improved normalization method, v2 regularization


# 3. other functions
1. visualization
2. computing cell cycle phase scores based on marker genes
3. DE testing
4. convert formats for different tools
5. parallelization to speed up

# 4. wrappers
fastMNN
Harmony
LIGER
Monocle3
Velocity


# STD workflow

library(dplyr)
library(Seurat)
library(patchwork)

############# generate object and QC
#############
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "../data/pbmc3k/filtered_gene_bc_matrices/hg19/") # genes are the first column; samples are column names
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# for mouse, use "^Mt-"
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Show QC metrics, e.g., metadata for the first 5 cells
head(pbmc@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, can be used for any columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# filteration based on QC metrcis using subset function
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

################# normalize, scale, and regress out
# normalization and log transformation: employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, 
# multiplies this by a scale factor (10,000 by default), and log-transforms the result.
# pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000) is the default and has been replaced by below 
pbmc <- NormalizeData(pbmc)

# finding HVG, e.g., feature selection
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# scaling the data using a linear transformation before dimensional redution: shift mean expression as 0, scale to variance across cells is 1 
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- ScaleData(pbmc) # the default only uses only 2000 HVG. In this case, the resulting DoHeatmap() will not cover all genes, but only the Scaled genes

# to remove/regress out unwated sources of variation, e.g, cell cycle stage, or mitochondrial contamination 
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt") # v3
# However, particularly for advanced users who would like to use this functionality, we strongly recommend the use of our new normalization workflow, SCTransform() with a vars.to.regress paramete
%%%%%%%%%%%%%%

######################### dimensional reduction and visualize it
# DEFAULT DIMENISONALITY IS 30, range can be 10 to 50
# linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5) # check the top and bottom 5 genes/features in the first dimension (PC1 TO PC5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca") # check the loading/weight of most contributing features (~ 30 genes) on PC1 and PC2
DimPlot(pbmc, reduction = "pca") # pca plot
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE) # heatmap for PC1, showing the clusters of 500 cells according to the exp of the most contributing features
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE) # heatmap for PC1 to PC15

# Determine the ‘dimensionality’ of the dataset using elbow plot and jackStraw plot
# NOTE: jackStraw method can take a long time for big datasets, comment out for expediency
#pbmc <- JackStraw(pbmc, num.replicate = 100)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
#JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)

######################## cluster
# first construct a KNN graph based on the euclidean distance in PCA space, 
# and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
#  Louvain algorithm (default) or SLM to iteratively group cells together,  
# resolution 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells. Optimal resolution often increases for larger datasets
pbmc <- FindClusters(pbmc, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

####################### Run non-linear dimensional reduction (UMAP/tSNE)
# As input to the UMAP and tSNE, we suggest using the same PCs as input to the clustering analysis
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages = 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(pbmc, reduction = "umap")


####################### Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 2
# The min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells
# The thresh.test argument requires a feature to be differentially expressed (on average) by some amount between the two groups, You can set both of these to 0, but with a dramatic increase in time
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive   ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC) # n = 2 means only get the value from the first 2 PCs

# the ROC test returns the ‘classification power’ for any individual marker (ranging from 0 - random, to 1 - perfect
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

####################### save
saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")


###################### visualization
# VlnPlot() show expression of features accross clusters
VlnPlot(pbmc, features = c("MS4A1", "CD79A")) 
# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

# FeaturePlot() visualizes feature expression on a tSNE or PCA plot
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
    "CD8A"))

# DoHeatmap() generates an expression heatmap for given cells and features
pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# RidgePlot(), CellScatter(), and DotPlot()


####################### Assigning cell type identity to clusters when we already know the corresponding cell type for each cluster

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


####################### save
saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")