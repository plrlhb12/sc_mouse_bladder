# load library
library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratObject)
library(ggplot2)
#options(Seurat.object.assay.version = "v5") # only need when annotate for ATAC, otherwise v4 is enough


######## haven't tested it yet

# Get the input argument so can run this script in terminal
# we need 3 input, 
# 1st is the input folder name or file name, 
# 2nd is the projectName, but the 2nd if optional but prefered in this script
# Rscript my_script.R "argument1" "argument2"

args <- commandArgs(trailingOnly = TRUE)

# check whether there are two input arguments
# the project_name will become the orig.ident
if (length(args) == 2) {
  project_name = args[2]
  cat("The input project name is:", args[1], "\n")
  cat("The input data is:", argsc, "\n")
} else {
  cat("The input data is:", args[1], "\n")
  cat("no project name was provided")
  project_name <- NULL
}

dataName <- args[1]
# Get information about the input argument whether it is a dir or a file
info <- file.info(args[1])

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

cat("sucessful")
cat(sessionInfo())

lung72h <- CreateSeuratObject(counts = lung72h.data, project = project_name, min.cells = 3, min.features = 200)
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
print("default Indents are:", Idents(lung72h))
Idents(lung72h) <- "seurat_clusters"
lung72h.markers <- FindAllMarkers(lung72h, min.pct = 0.25, logfc.threshold = 0.25) # only.pos = TRUE
write.csv(lung72h.markers, file = paste0(project_name, "_cluster_markers.csv"), row.names = FALSE)
saveRDS(lung72h, file = paste0(project_name, "std.rds"))

VlnPlot(lung72h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)# features = "percent.mt", split.by = "groups"
ggsave(paste0(Project, "_QC.pdf", plot2)# it autoclose the object, no object avail in console
    
plot1 <- FeatureScatter(lung72h, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(lung72h, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(paste0(Project, "_QC2.pdf"), (plot1+plot2))

plot3 <- ElbowPlot(lung72h)
ggsave(paste0(Project, "_elbow.pdf"), plot3)

plot4 <- DimPlot(lung72h, reduction = "pca") + NoLegend()# group.by = orig.ident or new annotation, split.by = conditions/groups in metadata
plot5 <- DimPlot(lung72h)
ggsave(paste0(Project, "_umap_pca.pdf"), (plot4 + plot5))
       
sink(paste0(project_name, "_log.txt")