# Data Integration: integrated analysis of single-cell datasets generated across different conditions, technologies, or species

############ integrate datasets in order to identfy and compare shared cell types across experiments
# rationale: identifying cell populations that are present across multiple datasets can be problematic under standard workflows

# methods: 
# first identify cross-dataset pairs of cells that are in a matched biological state (‘anchors’), 
# use anchors to correct for technical differences between datasets (i.e. batch effect correction), and to perform comparative scRNA-seq analysis of across experimental conditions.

list.dirs(path = "Y:/IRF/Cellranger/", recursive=FALSE)

# goals:
library(Seurat)
library(SeuratData)
library(patchwork)
InstallData("ifnb") # install dataset of human PBMC in resting or interferon-stiumulated states
LoadData("ifnb") # load dataset

########################## generate a list of datasets for integrating
# split the dataset into **a list** of two seurat objects (stim and CTRL)
# for real datasets, just create a list of them
ifnb.list <- SplitObject(ifnb, split.by = "stim")

########################## normalize and identify variable features for each dataset inside the dataset list independently
# option 1: std normalizatin
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration, this function was skipped in the mapp and annotate tutorial
features <- SelectIntegrationFeatures(object.list = ifnb.list) 
# for any integration in Seurat, need to use FindIntegrationAnchor and IntegrateData
# identify anchors, e.g., cell types present across datasets based on HVG
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
# creates an 'integrated' data assay based on identified anchors, here some batch correction or aligning was performed
# note that the original unmodified or uncorrected data still resides in the 'RNA' assay
immune.combined <- IntegrateData(anchorset = immune.anchors)

# option 2: improved normalization SCTransform
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
    anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

###############################
# run a std sc analysis on the integrated dataset on all cells!
DefaultAssay(immune.combined) <- "integrated" # swith default assay to the integrated dataset
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim") # umap, color by condtion
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE) # umap, color by clusters
p1 + p2
# or visualize condtions side by site on the sampe graph
DimPlot(immune.combined, reduction = "umap", split.by = "stim") # umap, color by clusters and split by conditions

# identify canonical cell type marker genes that are conserved across conditions
DefaultAssay(immune.combined) <- "RNA" # switch back to orignal data
nk.markers <- FindConservedMarkers(immune.combined, ident.1 = 6, grouping.var = "stim", verbose = FALSE) # compare cluster 6 to all others
head(nk.markers)
# plot expression of genes of interest on umap
FeaturePlot(immune.combined, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A","CCL2", "PPBP"), min.cutoff = "q9")
# manual change cell types according to pre-knowledge base
immune.combined <- RenameIdents(immune.combined, 
	`0` = "CD14 Mono", 
	`1` = "CD4 Naive T", 
	`2` = "CD4 Memory T",
    `3` = "CD16 Mono", 
    `4` = "B", 
    `5` = "CD8 T", 
    `6` = "NK", 
    `7` = "T activated", 
    `8` = "DC", 
    `9` = "B Activated",
    `10` = "Mk", 
    `11` = "pDC", 
    `12` = "Eryth", 
    `13` = "Mono/Mk Doublets", 
    `14` = "HSPC")
DimPlot(immune.combined, label = TRUE)

# plot dot plot which showing both expression level and percentage
# need to change characters of Idents(obj) to factors, we can define the sequences of the cell types showing in the plot
# by using levels = c()
Idents(immune.combined) <- factor(Idents(immune.combined), levels = c("HSPC", "Mono/Mk Doublets",
    "pDC", "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono", "B Activated", "B", "CD8 T", "NK", "T activated",
    "CD4 Naive T", "CD4 Memory T"))
# define the genes to plot, picking 2-3 strong maker genes for all the clusters
markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5",
    "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",
    "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", "PRSS57")
# dot plot
DotPlot(immune.combined, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "stim") +
    RotatedAxis()


######### Identify and visual differential expressed genes across conditions
# plot the average expression of both the stimulated and control cells and look for genes that are visual outliers on a scatter plot
# only examine two cells types of interest for comparisons in this example
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
# extract CD Navie T cells
t.cells <- subset(immune.combined, idents = "CD4 Naive T")
Idents(t.cells) <- "stim" # assing the info of stim column as its indentity
avg.t.cells <- as.data.frame(log1p(AverageExpression(t.cells, verbose = FALSE)$RNA))
avg.t.cells$gene <- rownames(avg.t.cells)
# extract CD14 Mono cells
cd14.mono <- subset(immune.combined, idents = "CD14 Mono")
Idents(cd14.mono) <- "stim" # assing the info of stim column as its indentity
avg.cd14.mono <- as.data.frame(log1p(AverageExpression(cd14.mono, verbose = FALSE)$RNA))
avg.cd14.mono$gene <- rownames(avg.cd14.mono)
# define the genes to be labeled on the scatter plot
genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")
# scatter plot on CD4 Navie T cells
p1 <- ggplot(avg.t.cells, aes(CTRL, STIM)) + geom_point() + ggtitle("CD4 Naive T Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
# scatter plot on CD14 Monocytes
p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
p1 + p2

# Because we are confident in having identified common cell types across condition, 
# we can ask what genes change in different conditions for cells of the same type
# what need to do its generate new values of Idents(obj) based on contrl and stimultaion and perform DE as std workflow
# add one column of "celltype.stim" info in the meta.data slot, which containing cellTypeIdentity and stim info
immune.combined$celltype.stim <- paste(Idents(immune.combined), immune.combined$stim, sep = "_")
# store the orignal identiy as another column in meta.data
immune.combined$celltype <- Idents(immune.combined)
# switch the identify column's info
Idents(immune.combined) <- "celltype.stim"
# perform DE on conditions instead of on clusters
b.interferon.response <- FindMarkers(immune.combined, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)
head(b.interferon.response, n = 15)
# plot gene expression between conditions on either umap scatter or violn graph
FeaturePlot(immune.combined, features = c("CD3D", "GNLY", "IFI6"), split.by = "stim", max.cutoff = 3,
    cols = c("grey", "red"))
plots <- VlnPlot(immune.combined, features = c("LYZ", "ISG15", "CXCL10"), split.by = "stim", group.by = "celltype",
    pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

################ Tips for integrating very large dataset
# https://satijalab.org/seurat/articles/integration_large_datasets.html 
# updated on 20240124: this page was not avail anymore. I executed my succesfully on Sep2023; Seurat updated its vignettes on Dec2023. 
# When I troubleshoot my integration between allCTRL AND allviuse this page didn't exist anymore
In general, we observe strikingly similar results between the standard workflow and the one demonstrated here, with substantial reduction in compute time and memory. 
However, if the datasets are highly divergent (for example, cross-modality mapping or cross-species mapping), where only a small subset of features can be used to facilitate integration, and you may observe superior results using CCA.
# two options
# Option 1: RPCA: reciprocal PCA
# we use reciprocal PCA (RPCA) instead of CCA, to identify an effective space in which to find anchors in in FindIntegrationAnchors()
# Option 2: reference-based integration
# instead of all pairing, specifying one or mores as reference for integration while the remainder as query for integration
library(Seurat)
bm280k.data <- Read10X_h5("../data/ica_bone_marrow_h5.h5")
bm280k <- CreateSeuratObject(counts = bm280k.data, min.cells = 100, min.features = 500)
bm280k.list <- SplitObject(bm280k, split.by = "orig.ident")
bm280k.list <- lapply(X = bm280k.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})

#  SelectIntegrationFeatures is required for running the alternative reciprocal PCA workflow.
features <- SelectIntegrationFeatures(object.list = bm280k.list)
bm280k.list <- lapply(X = bm280k.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = bm280k.list, reference = c(1, 2), reduction = "rpca",
    dims = 1:50) # only chose dataset 1 and 2 to use in a reference-based wf
bm280k.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
bm280k.integrated <- ScaleData(bm280k.integrated, verbose = FALSE)
bm280k.integrated <- RunPCA(bm280k.integrated, verbose = FALSE)
bm280k.integrated <- RunUMAP(bm280k.integrated, dims = 1:50)
DimPlot(bm280k.integrated, group.by = "orig.ident")


