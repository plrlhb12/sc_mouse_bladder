# see the functions of Seurat
# https://satijalab.org/seurat/reference/index.html

# quick cheatsheet
library(Seurat)
library(SeuratData)
library(dplyr)
library(sctransform)
library(patchwork)
library(Azimuth)
library(glmGamPoi) # for better performance of SCTtransform
library(BPCells) #
#library(scater) # for more options of visualization
#library(ggplot2)
#library(SeuratDisk) # for interoperationality
# options(Seurat.object.assay.version = "v5") # only optional

###########################################
######### my way to read data
###########################################
file_name <- "your_file.h5ad"

if (grepl("\\.h5ad$", file_name)) {
  print("The file name ends with '.h5ad'")
  obj <- RunAzimuth(file_name, reference="lungref")
} else if (grepl("\\.h5$", file_name)) {
  print("The file name ends with '.h5'")
  obj <- Read10X_h5(file_name)
  obj <- RunAzimuth(obj, reference="lungref")
} else if (file.info(file_name)$isdir) {
  print("The file name is a directory")
  obj <- Read10X(file_name)
  obj <- RunAzimuth(obj, reference="lungref")
} else {
  print("The file name does not match any condition")
}


###########################################
# read different input formats
###########################################

# input
Read10X() # for output dir from Cellranger
Read10X_h5() # for .h5 fle from Cellranger
readRDS() # for .rds saved an R object from Seurat
pbmc3k <- LoadData("pbmc3k") # from build-in dataset
hamsterlung <- readRDS("../../seurat_biowulf/FixedRNAPilot_seurat/lung72h.rds")
RunAzimuth(file.h5ad, reference="lungref")


###########################################
# std wf
# Part 1: std seurat wf without visualization
###########################################

Path = "seurat_biowulf/FixedRNAPilot_seurat"
Project = "Lung72h"
Counts <- Read10X(data.dir = "/data/pengl7/IRF/Cellranger/FixedRNAPilot9/outs/per_sample_outs/Lung72h_BC2/count/sample_filtered_feature_bc_matrix")

lung72h <- CreateSeuratObject(counts = Counts, project = Project, min.cells = 3, min.features = 200)
lung72h[["percent.mt"]] <- PercentageFeatureSet(lung72h, pattern = "^MT-") # change to ^Mt- for mouse, can skip this for 10x Flex plotform
VlnPlot(lung72h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(lung72h, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(lung72h, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
lung72h <- subset(lung72h, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
lung72h <- NormalizeData(lung72h)
lung72h <- FindVariableFeatures(lung72h, selection.method = "vst", nfeatures = 2000)
lung72h <- ScaleData(lung72, vars.to.regress = "percent.mito") # e.g., vars.to.regress = "percent.ERCC", ... etc
lung72h <- RunPCA(lung72h, features = VariableFeatures(object = lung72h))
ElbowPlot(lung72h) # to define the number of PCA for use
lung72h <- FindNeighbors(lung72h, dims = 1:10) # default dims = 1:10
lung72h <- RunUMAP(lung72h, dims = 1:10) # default dims = NULL
lung72h <- FindClusters(lung72h, resolution = 0.5)
lung72h <- RunAzimuth(lung72h, reference="lungref")
saveRDS(lung72h, file = paste0(project_name, "az.rds")
###########################################
# Part I: std seurat wf using a pipe function
# successfuly in biowulf R4.3 version in one command; 
# couldn't add plots more than 1 since the returned object will not be the Seurat Object anymore
# couldn't add FindMarker() and plot in pip, because the default seurat object changes after one of these functions
# thus add some functions outside the pipe function
# change pattern = "^MT-" to "^Mt-" for mouse
###########################################

Path = "seurat_biowulf/FixedRNAPilot_seurat"
Project = "lung72h30June"
Counts <- Read10X(data.dir = "/data/pengl7/IRF/Cellranger/FixedRNAPilot9/outs/per_sample_outs/Lung72h_BC2/count/sample_filtered_feature_bc_matrix")

lung72h <- CreateSeuratObject(Counts, project = Project, min.cells = 3, min.features = 200) %>%
    PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
    subset(subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) %>% ???? double check
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(vars.to.regress = "percent.mt") %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:10) %>%
    RunUMAP(dims = 1:10) %>%
    FindClusters(resolution = 0.5) %>%
    RunAzimuth(reference="lungref") %>%
    saveRDS(file = paste0(Project, "std_az.rds"))

##############################################
# Part II: find markers and visualize after pipe function
############################################## 
print("default Indents are:", Idents(lung72h))
Idents(lung72h) <- "seurat_clusters"
lung72h.markers <- FindAllMarkers(lung72h, min.pct = 0.25, logfc.threshold = 0.25) # only.pos = TRUE
write.csv(lung72h.markers, file = paste0(project_name, "_cluster_markers.csv"), row.names = FALSE)
saveRDS(lung72h, file = paste0(Project, ".rds"))

# features = "percent.mt", split.by = "groups"
VlnPlot(lung72h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# # it autoclose the object, no object avail in console
ggsave(paste0(Project, "_QC.pdf", plot2)
    
plot1 <- FeatureScatter(lung72h, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(lung72h, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(paste0(Project, "_QC2.pdf"), (plot1+plot2))

plot3 <- ElbowPlot(lung72h)
ggsave(paste0(Project, "_elbow.pdf"), plot3)

plot4 <- DimPlot(lung72h, reduction = "pca") + NoLegend()
# group.by = orig.ident or new annotation, split.by = conditions/groups/... in metadata
plot5 <- DimPlot(lung72h)
ggsave(paste0(Project, "_umap_pca.pdf"), (plot4 + plot5))



###########################################
# SCT version seurat wf
# Part 1:  without visualization
# SCTransform workflow, this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures().
# more options in impoving, e.g., using v2, and including glmGamPoi
# object <- SCTransform(object, vst.flavor = "v2") # using V2
##########################################
pbmc_data <- Read10X(data.dir = "../data/pbmc3k/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc_data, project = Project, min.cells = 3, min.features = 200)
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
# using glmGamPoi package which substantially improves the speed of the learning procedure
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("glmGamPoi")
pbmc <- SCTransform(pbmc, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
# These are now standard steps but using more default dims in the Seurat workflow for visualization and clustering
pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE)
lung72h <- RunAzimuth(lung72h, reference="lungref")
DimPlot(pbmc, label = TRUE) + NoLegend()
saveRDS(pbmc, file = paste0(Project, "sct_az.rds"))

## in a single command
# # successful in biowulf
Path = "seurat_biowulf/FixedRNAPilot_seurat"
Project = "lung72h_sct"
Counts <- Read10X(data.dir = "/data/pengl7/IRF/Cellranger/FixedRNAPilot9/outs/per_sample_outs/Lung72h_BC2/count/sample_filtered_feature_bc_matrix")

lung72h_sct <- CreateSeuratObject(Counts, project = Project, min.cells = 3, min.features = 200) %>%
    PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
    SCTransform(vars.to.regress = "percent.mt") %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>%
    RunUMAP(dims = 1:30) %>%
    FindClusters() %>%
    RunAzimuth(reference="lungref") %>%
    saveRDS(file = paste0(Project, "_az.rds"))


##############################################
# Part II: find markers and visualize after pipe function
############################################## 
print("default Indents are:", Idents(lung72h))
Idents(lung72h) <- "seurat_clusters" # change it according the clusters or cell types of interest
lung72h.markers <- FindAllMarkers(lung72h, min.pct = 0.25, logfc.threshold = 0.25) # only.pos = TRUE
write.csv(lung72h.markers, file = paste0(project_name, "_cluster_markers.csv"), row.names = FALSE)
saveRDS(lung72h, file = paste0(Project, ".rds"))

# features = "percent.mt", split.by = "groups"
VlnPlot(lung72h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# # it autoclose the object, no object avail in console
ggsave(paste0(Project, "_QC.pdf", plot2)
    
plot1 <- FeatureScatter(lung72h, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(lung72h, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(paste0(Project, "_QC2.pdf"), (plot1+plot2))

plot3 <- ElbowPlot(lung72h)
ggsave(paste0(Project, "_elbow.pdf"), plot3)

plot4 <- DimPlot(lung72h, reduction = "pca") + NoLegend()
# group.by = orig.ident or new annotation, split.by = conditions/groups/... in metadata
plot5 <- DimPlot(lung72h)
ggsave(paste0(Project, "_umap_pca.pdf"), (plot4 + plot5))


############################################################
# save
############################################################
saveRDS(lung72h, file = "lung72h.rds")
ggsave(filename=paste0(Project, "_pca.pdf"), (plot1+plot2), width = 6, height = 6, dpi = 300)
ggsave("dimplot_lung72h.pdf", width = 14, height = 7, units = "in", DimPlot(lung72, pt.size = 3, label = TRUE))


############################################################
# annotaiton using Azimuth
# Notice that the default Assay becomes "refRNA", however, the counts layer inside refRNA is the same as in "RNA"
############################################################

AvailableData()
available_data[grep("Azimuth", available_data[, 3]), 1:3]
lung72 <- readRDS("lung72h.rds") # if starting from beginning
lung72 <- RunAzimuth(lung72, reference = "lungref")
DefaultAssay(lung72h)  #"refAssay"
Assays(lung72h)
DefaultAssay(lung72h) <- "RNA" #if need to switch

> colnames(lung72h@meta.data)
[1] "orig.ident"                      
[2] "nCount_RNA"                      
[3] "nFeature_RNA"                    
[4] "percent.mt"                      
[5] "RNA_snn_res.0.5"                 
[6] "seurat_clusters"                 
[7] "nCount_refAssay"                 
[8] "nFeature_refAssay"               
[9] "predicted.ann_level_1.score"     
[10] "predicted.ann_level_1"           
[11] "predicted.ann_level_2.score"     
[12] "predicted.ann_level_2"           
[13] "predicted.ann_level_3.score"     
[14] "predicted.ann_level_3"           
[15] "predicted.ann_level_4.score"     
[16] "predicted.ann_level_4"           
[17] "predicted.ann_level_5.score"     
[18] "predicted.ann_level_5"           
[19] "predicted.ann_finest_level.score"
[20] "predicted.ann_finest_level"      
[21] "mapping.score"  

#######################################################################################################################
# find markers; the default groups used for DE will be the Idents(), need to chagne its value fo findMarkers() to work
# check the Indents() whether the grouping is of interet; otherwise may error out
# need to assign to a dataframe oject then write.csv() otherwise it will not saved in the Seurat object
#######################################################################################################################
# pbmc.markers %>% group_by(cluster)
FindAllMarkers(lung72h, min.pct = 0.25, logfc.threshold = 0.25)

# possible error: Cell group 2 is empty - no cells with identity class
head(lung72h@meta.data, 3)
FindAllMarkers(lung72h, min.pct = 0.25, logfc.threshold = 0.25)
colnames(lung72h@meta.data)
levels(lung72h$seurat_clusters)
table(lung72h$seurat_clusters)
table(as.factor(lung72h$seurat_clusters))
table(as.factor(lung72h$predicted.ann_level_1))
head(Idents(lung72h),10)  # the original indents are cell barcodes

DefaultAssay(lung72h) <- "RNA" # doesn't need to swith, the count layer are the same inside "refRNA" and "RNA" ???
allMarkers <- FindAllMarkers(min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(allMarkers, file = paste0(Project, "_markers.csv"), row.names = FALSE)
Idents(lung72h) <- "seurat_clusters"
markers <- FindAllMarkers(lung72h, min.pct = 0.25, logfc.threshold = 0.25)
head(markers)

pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC) # n = 2 means only get the value from the first 2 PCs


################# after one click generation of object
names(lung72) # show assay names
Idents(lung72) # check the default name of each cells, need to change to desired such as cluster numbers of annotated cell type names
colnames(lung72@meta.data) # show column names of metadata dataframe
levels(lung72$orig.ident) # check the numbers of the levels