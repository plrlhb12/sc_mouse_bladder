# load library
library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratObject)
library(ggplot2)
library(SCTransform)
library(Azimuth)
library(glmGamPoi)
#options(Seurat.object.assay.version = "v5") # only need when annotate for ATAC, otherwise v4 is enough

dataName <- "....."
project_name <- "...."

# Get information about the input argument whether it is a dir or a file
info <- file.info(dataName)

if (info$isdir){
  cat(dataName, "is a directory. \n")
  counts <- Read10X(data.dir = dataName)
} else if (endsWith(dataName, ".h5")) {
  counts <- Read10X_h5(dataName)
} else if (endsWith(dataName, ".h5ad")) {
  counts <- ReadRDS(dataName)    
} else {
  stop(paste0("Error: ", dataName, " does not exist."))
}

################################## std
lung72h <- CreateSeuratObject(counts = counts, project = project_name, min.cells = 3, min.features = 200)
lung72h[["percent.mt"]] <- PercentageFeatureSet(lung72h, pattern = "^Mt-")
VlnPlot(lung72h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(lung72h, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(lung72h, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
lung72h <- subset(lung72h, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
lung72h <- NormalizeData(lung72h)
lung72h <- FindVariableFeatures(lung72h, selection.method = "vst", nfeatures = 2000)
lung72h <- ScaleData(lung72) # e.g., vars.to.regress = "percent.mito" or "percent.ERCC", ... etc
lung72h <- RunPCA(lung72h, features = VariableFeatures(object = lung72h))
ElbowPlot(lung72h) # to define the number of PCA for use
lung72h <- FindNeighbors(lung72h, dims = 1:10) # default dims = 1:10
lung72h <- RunUMAP(lung72h, dims = 1:10) # default dims = NULL
lung72h <- FindClusters(lung72h, resolution = 0.5)
lung72h <- RunAzimuth(lung72h, reference="lungref")
print("default Indents are:", Idents(lung72h))
Idents(lung72h) <- "seurat_clusters"
lung72h.markers <- FindAllMarkers(lung72h, min.pct = 0.25, logfc.threshold = 0.25) # only.pos = TRUE
write.csv(lung72h.markers, file = paste0(project_name, "_cluster_markers.csv"), row.names = FALSE)
saveRDS(lung72h, file = paste0(project_name, "std_az.rds"))


########################################## sct
pbmc_data <- Read10X(data.dir = "../data/pbmc3k/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc_data, project = Project, min.cells = 3, min.features = 200)
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
# using glmGamPoi package which substantially improves the speed of the learning procedure
#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("glmGamPoi")
pbmc <- SCTransform(pbmc, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
# These are now standard steps but using more default dims in the Seurat workflow for visualization and clustering
pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE)
lung72h <- RunAzimuth(lung72h, reference="lungref")
DimPlot(pbmc, label = TRUE) + NoLegend()
print("default Indents are:", Idents(lung72h))
Idents(lung72h) <- "seurat_clusters"
lung72h.markers <- FindAllMarkers(lung72h, min.pct = 0.25, logfc.threshold = 0.25) # only.pos = TRUE
write.csv(lung72h.markers, file = paste0(project_name, "_cluster_markers.csv"), row.names = FALSE)
saveRDS(pbmc, file = paste0(project_name, "sct_az.rds"))
