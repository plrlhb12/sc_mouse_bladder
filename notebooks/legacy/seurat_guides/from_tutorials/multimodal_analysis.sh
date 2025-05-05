# https://satijalab.org/seurat/articles/multimodal_vignette.html
# this example clusters a CITE-seq dataset on the basis of the measured cellular transcriptomes, and subsequently discover cell surface proteins that are enriched in each cluster
# sample dataset:  8,617 cord blood mononuclear cells (CBMCs), where transcriptomic measurements are paired with abundance estimates for 11 surface proteins quantified with DNA-barcoded antibodies
# more advanced: the application of our Weighted Nearest Neighbors (WNN) approach that enables simultaneous clustering of cells based on a ** weighted combination of both modalities **


# cluster CITE on transcriptiomes

library(Seurat)
library(ggplot2)
library(patchwork)

# load 2 module's matrix seperately and make sure they have the same cells
# Load in the RNA UMI matrix, and ADT UMI matrix seperately
cbmc.rna <- as.sparse(read.csv(file = "../data/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", sep = ",",
    header = TRUE, row.names = 1))
# Since this dataset has mouse and human genes, make life a bit easier going forward, we're going to discard all but the top 100 most highly expressed mouse genes as negative control, and remove the 'HUMAN_' from the CITE-seq prefix
cbmc.rna <- CollapseSpeciesExpressionMatrix(cbmc.rna)

cbmc.adt <- as.sparse(read.csv(file = "../data/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz", sep = ",",
    header = TRUE, row.names = 1))

# Note that since measurements were made in the same cells, the two matrices have identical column names
all.equal(colnames(cbmc.rna), colnames(cbmc.adt))


# set a Seurat object containing two types of assay, one for RNA, the other one for antibody-derived tags (ADT)
cbmc <- CreateSeuratObject(counts = cbmc.rna)
Assays(cbmc) # first generate seurat object using sc data
adt_assay <- CreateAssayObject(counts = cbmc.adt) # genearate an Assay object and assign the above Seurate object
cbmc[["ADT"]] <- adt_assay
rownames(cbmc[["ADT"]]) # check the features of CITE, here are 11 antibodies's names
DefaultAssay(cbmc)
DefaultAssay(cbmc) <- "ADT" # change default assay when necessary
DefaultAssay(cbmc)

# Cluster cells on the basis of their scRNA-seq profiles
# perform visualization and clustering steps
DefaultAssay(cbmc) <- "RNA"
DefaultAssay(cbmc)
cbmc <- NormalizeData(cbmc)
cbmc <- FindVariableFeatures(cbmc)
cbmc <- ScaleData(cbmc)
cbmc <- RunPCA(cbmc, verbose = FALSE)
cbmc <- FindNeighbors(cbmc, dims = 1:30)
cbmc <- FindClusters(cbmc, resolution = 0.8, verbose = FALSE)
cbmc <- RunUMAP(cbmc, dims = 1:30)
DimPlot(cbmc, label = TRUE) # show clusters on umap

############################
# Visualize multiple modalities side-by-side; 
# Method 1: switch default assay back and forth to perform functions on the right assay of interest
DefaultAssay(cbmc) <- "ADT"
cbmc <- NormalizeData(cbmc, normalization.method = "CLR", margin = 2) # normalize ADT data
# Note that the following command is an alternative but returns the same result: cbmc <- NormalizeData(cbmc, normalization.method = "CLR", margin = 2, assay = "ADT")
p1 <- FeaturePlot(cbmc, "CD19", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
DefaultAssay(cbmc) <- "RNA"
p2 <- FeaturePlot(cbmc, "CD19") + ggtitle("CD19 RNA")
p1 | p2 # place plots side-by-side


# Method 2: define the key for the assay of interest inside functions
Key(cbmc[["RNA"]]) ## [1] "rna_"
Key(cbmc[["ADT"]]) ## [1] "adt_"
p1 <- FeaturePlot(cbmc, "adt_CD19", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein") # using adt_ to define the assay is ADT
p2 <- FeaturePlot(cbmc, "rna_CD19") + ggtitle("CD19 RNA")
p1 | p2

########################### visual expression and identify cluster markers
# leverage our paired CITE-seq measurements to help annotate clusters derived from scRNA-seq, and to identify both protein and RNA markers
# as we know that CD19 is a B cell marker, we can identify cluster 6 as expressing CD19 on the surface
VlnPlot(cbmc, "adt_CD19") # plot violin on the assay of ADT using CD19 expression, it shows that cluster 6 expression high level of CD19, so cluster 6 is B cell

# identify alternative protein and RNA markers for this cluster through differential expression
adt_markers <- FindMarkers(cbmc, ident.1 = 6, assay = "ADT")
rna_markers <- FindMarkers(cbmc, ident.1 = 6, assay = "RNA")
head(adt_markers)
head(rna_markers)

########################### addtional visualization
# Draw ADT scatter plots (like biaxial plots for FACS). Note that you can even 'gate' cells if desired by using HoverLocator and FeatureLocator
FeatureScatter(cbmc, feature1 = "adt_CD19", feature2 = "adt_CD3")
# view relationship between protein and RNA
FeatureScatter(cbmc, feature1 = "adt_CD3", feature2 = "rna_CD3E")
FeatureScatter(cbmc, feature1 = "adt_CD4", feature2 = "adt_CD8")

########################### Loading data from 10X multi-modal experiments freely available from 10X Genomics
pbmc10k.data <- Read10X(data.dir = "../data/pbmc10k/filtered_feature_bc_matrix/") # read count dir
# do some changes on pbmc10k.data[["Antibody Capture"]]'s rownames, i.e., antibody names
rownames(x = pbmc10k.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "", 
x = rownames(x = pbmc10k.data[["Antibody Capture"]]))
pbmc10k <- CreateSeuratObject(counts = pbmc10k.data[["Gene Expression"]], min.cells = 3, min.features = 200)
pbmc10k <- NormalizeData(pbmc10k)
pbmc10k[["ADT"]] <- CreateAssayObject(pbmc10k.data[["Antibody Capture"]][, colnames(x = pbmc10k)])
pbmc10k <- NormalizeData(pbmc10k, assay = "ADT", normalization.method = "CLR")
plot1 <- FeatureScatter(pbmc10k, feature1 = "adt_CD19", feature2 = "adt_CD3", pt.size = 1)
plot2 <- FeatureScatter(pbmc10k, feature1 = "adt_CD4", feature2 = "adt_CD8a", pt.size = 1)
plot3 <- FeatureScatter(pbmc10k, feature1 = "adt_CD3", feature2 = "CD3E", pt.size = 1)
(plot1 + plot2 + plot3) & NoLegend()

# Additional functionality for multimodal data in Seurat
# Defining cellular identity from multimodal data using WNN analysis in Seurat v4 vignette
# Mapping scRNA-seq data onto CITE-seq references [vignette]
# Introduction to the analysis of spatial transcriptomics analysis [vignette]
# Analysis of 10x multiome (paired scRNA-seq + ATAC) using WNN analysis [vignette]
# Signac: Analysis, interpretation, and exploration of single-cell chromatin datasets [package]
# Mixscape: an analytical toolkit for pooled single-cell genetic screens [vignette]