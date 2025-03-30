# most important things
DefaultyAssay(obj) <- "RNA" # for DE and AVG expr
DefaultyAssay(obj) <- "integrated" # for clustering and plotting
Indent(obj) <- "cell.type" # need if to define idents inside the FindMarker function


# input, e.g., rds, h5, an count folder
readRDS()
saveRDS()

# file, folder, and paths
setwd("Y:/IRF/Cellranger/seurat_biowulf/FixedRNAPilot_seurat")
getwd()
list.dirs(recursive = FALSE)
list.dirs(path = "../pool1_nextseq/nextseq2/per_sample_outs", recursive=FALSE)
list.files("../pool1_nextseq/nextseq2/per_sample_outs/065E-1B-Epsilon-D10/count")
list.files(path=PATH, pattern = "Ctrl.*rds$") # only list files containing Ctrl anad end with rds
list.files(path=PATH, pattern = "\\.rds$") # only list files end with .rds
run.name <- "1B"
list.files(path = "h5files", pattern = run.name)
# Create the new folder
dir.create(run.name)
# Check if the folder was created
if (file.info(run.name)$isdir) {
  cat("Folder created successfully.\n")
} else {
  cat("Failed to create folder.\n")
}
# get folder name
dirname(file)
# make a path: foldername + filename
new_path <- file.path(dirname(file), new_name)
# rename a file or folder
file.rename(old_name, new_path)


# generate a list of marker genes
# gene <- "Siglecf,Marco,C1qb,Ccr2,Ccr5,Arg1,Treml4,Foxj1,S100a8,Cxcr2,Camp,Flt3,H2-Ab1hi,Irf8lo,Tcf4lo,Flt3,H2-Ab1hi,Irf8hi,Tcf4hi,d3e,Cd4orCd8a,Gzma,Nkg7,Cd79b,Ms4a1"
elements <- strsplit(gene, ",")[[1]]
gene_list <- c(elements)
gene_list
DefaultAssay(immune) <- "RNA"
FeaturePlot(immune, features = gene_list)
ggsave("4J-Omicron-D3/20230818/immune_makers_in_subset.pdf")

# remove A sublist from a list, i..e, subset a list
# using boolean list to subset
whole <- list.files(path = "../20230823-individual/individuals-2/", pattern = ".*\\.rds$") # get A list containg file names
viral <- whole %in% ctrl.list #get a boolean list which False/True in B ctrl list
virus <- whole[!viral] # subset the A big list using the boolean list
# using grep to subset
d3.list <- virus[grep(pattern = ".*D3.rds")]


# string operations using pattern or regular expression 
# The double backslash (\\) before the dot (.) is used to escape the dot because 
# in regular expressions, a dot has a special meaning (matches any character)
batches <- sub("^065E-4J-(.*?)\\.h5$", "\\1", h5_file_list)
batches <- gsub("\\.h5", "", h5_file_list)
batches <- gsub("^065E-|\\.h5$", "", h5_file_list)
obj$orig.ident <- gsub("oldname", "newname", obj$orig.ident)
sub.obj@meta.data[["replicates"]] <- sub("^...", "", sub.obj@meta.data$orig.ident) # replace the first 3 charactes with an empty string
sub.obj@meta.data[["strains"]] <- sub("-.*", "", sub.obj@meta.data$replicates) # replace 
obj@meta.data[["strains.fine.level"]] <- sub("....$", "", obj$strains.fine.level)
new_i <- gsub(pattern = "^2023-08-28-|^2023-08-29-|\\.rds$", "", i) # replace "2020-0829-" or "20230829-" or ".rds" with ""


SessionInfo()

library(Seurat)
library(dplyr)
library(patchwork)
library(SeuratData)
library(sctransform)
library(Azimuth)
library(ggplot2)
library(hdf5r)

# calculate the frequence of the categories in a column, e.g., calculate cell number per type
table(obj@meta.data$days)

sink("allvirus_sars_cellType_day.txt")
frequency_table <- table(meta$cell.type, meta$days, meta$sars.pos)
print(frequency_table)
sink()

# if want to save the std out to an execel file and covert a long form to a wide form
a.table <- table(sub.obj$predicted.celltype, sub.obj$days)
df <- data.frame(a.table)
df <- pivot_wider(df, id_cols = Var1, names_from = Var2, values_from = Freq)
write.csv(df, file = paste0(run.name, "/", run.name, "_cell_numbers_per_type.csv"), row.names = FALSE)

levels(lung72h$seurat_clusters)
levels(factor(sub.obj@meta.data$orig.ident))
table(as.factor(lung72h$seurat_clusters))
table(as.factor(lung72h$predicted.ann_level_1))
head(Idents(lung72h),10)  # the original indents are cell barcodes
Idents(lung72h) <- "seurat_clusters"
markers <- FindAllMarkers(lung72h, min.pct = 0.25, logfc.threshold = 0.25)
head(lung72h@meta.data, 3)
DefaultAssay(lung72h)  #"refAssay"
Assays(lung72h)
DefaultAssay(lung72h) <- "RNA"
FindAllMarkers(lung72h, min.pct = 0.25, logfc.threshold = 0.25)
colnames(lung72h@meta.data)

# check the expression of SARS
"M" %in% rownames(combined)
"ORF1ab" %in% rownames(d3[["RNA"]])

sars1 <- c("S", "M", "N", "ORF1ab", "ORF8")
genes_present <- is.element(sars1, genes$V2)

# check the original feature table
genes <- read.delim("../features.tsv", header = FALSE)
dim(genes)


########## subset
###################### subset a dataframe using boolean condition
significant_genes <- de_results[de_results$p_val_adj < 0.05 & abs(de_results$log2fc) > 1, ]


Idents(combined) <- combined$orig.ident
d3 <- subset(combined, idents = "D3")
Idents(d3) <- d3$seurat_clusters

lung72@assays[["refAssay"]]
names(lung72)
colnames(lung72@meta.data)
lung72@assays[["refAssay"]][1:5, 1:2] 
lung72@assays[["ref.umap"]]
lung72@assays[["ref.umap"]]

# subset immune cells and re-clustering
Idents(d3) <- d3$predicted.celltype
types.remove <- c("AT1", "AT2", "Ciliated", "Endothelial", "Fibroblasts", "Unclear", "SmoothMuscle", "Myofibroblast")
immune <- subset(d3, idents = types.remove, invert = TRUE)

# since the color are so simlliar, show the types separately
for (type in levels(factor(d3$predicted.celltype))){
  DimPlot(d3, group.by = "predicted.celltype", cells = d3$predicted.celltype == type) +xlim(-15,15) +ylim(-15,15)
  ggsave(paste0("4J-Omicron-D3/20230818/", type, "_umap.pdf"))
}


########## reference data
AvailableData()
InstallData("lungref")
data("lungref")
RemoveData('lungref')
available_data[grep("Azimuth", available_data[, 3]), 1:3]
lung_ref <- LoadData("lungref", "azimuth")


##################################

# raster parameter problem when too many points (exceeding 100,000) in a large dataset in plotting
# pt.size = 1 Sset raster set itself automatically based on object size
# using `raster=FALSE` has better visual resolutions, the output pdf is very big which may take a minute or longer to open it

DimPlot(sub.obj, group.by = "predicted.celltype", pt.size = 1)

# a whole pipeline for single dataset
# ggave
p3 <- DimPlot(sub.obj, group.by = "predicted.celltype", label = TRUE, pt.size = 1, label.size = 2) + NoLegend() #pt.size = 1 Sset raster set itself automatically based on object size.
p1 + p2 + p3
ggsave(filename = paste0(run.name, "/", "/integrated_virus.pdf"), plot = p3, width = 16)

lung72 <- readRDS("seurat_biowulf/FixedRNAPilot_seurat/lung72h_new.rds")
lung72
VlnPlot(lung72, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(lung72, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(lung72, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave(filename=paste0("qc", "_.pdf"), (plot1+plot2), width = 6, height = 6, dpi = 300)

ElbowPlot(combined)

FeatureScatter(combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident") # MUST change group.by if obj has been clustered

# DoHeatmap need scaled data
#d3 <- ScaleData(d3)
DefaultAssay(d3) <- "SCT" # SCT can be used for visualization, however be cautions if there are seveal dataset together with big difference of sequencing depth
DoHeatmap(d3, features = top10_cluster_markers$gene, size = 5)+NoLegend() + theme(text = element_text(size = 8))
ggsave(filename = "4J-Omicron-D3/20230818/4J-Omicron-D3-heatmap_top10_18clusters.pdf", height = 14, width = 12)

# sct pipeline
# remember to change the values of subset(subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# for ngeative bionomial regress only using 2000 genes 5000 cells
project = "aBatchName"
lung72h.data <- Read10X(data.dir = "/data/pengl7/IRF/Cellranger/FixedRNAPilot9/outs/per_sample_outs/Lung72h_BC2/count/sample_filtered_feature_bc_matrix")
lung72h_sct <- CreateSeuratObject(lung72h.data, project = project) %>%
  PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
  #subset(subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) %>% I didn't add this line when exploring the pilot data
  SCTransform(vars.to.regress = "percent.mt") %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  RunUMAP(dims = 1:30) %>%
  FindClusters() %>%
  RunAzimuth(reference="lungref") %>%
  saveRDS(file = "lung72h_sct.rds")

allMarkers <- FindAllMarkers(lung72h-sct, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(allMarkers, file = paste0(Project, "_markers.csv"), row.names = FALSE)

# sd workflow
# remember to change the values of subset(subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
project = "aBatchName"
Counts <- Read10X(data.dir = "/data/pengl7/IRF/Cellranger/FixedRNAPilot9/outs/per_sample_outs/Lung72h_BC2/count/sample_filtered_feature_bc_matrix")
lung72h-test <- CreateSeuratObject(Counts, project = project) %>%
  PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
  subset(subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5) %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(vars.to.regress = "percent.mt") %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:10) %>%
  RunUMAP(dims = 1:10) %>%
  FindClusters(resolution = 0.5) %>%
  RunAzimuth(reference="lungref") %>%
  saveRDS(file = paste0(Project, ".rds"))

# export some makers
DefaultAssay(d3) <- "RNA"
cluster21_makers <- FindAllMarkers(d3, only.pos = TRUE)
write.csv(cluster21_makers, file = "4J-Omicron-D3/20230818/cluster18_marker_genes_4J-Omicron-D3.csv", row.names = FALSE)


######################################
# integrate several 
######################################### 
list.dirs(path = "../pool1_nextseq/nextseq2/per_sample_outs", recursive=FALSE)
list.files("../pool1_nextseq/nextseq2/per_sample_outs/065E-1B-Epsilon-D10/count")
h5_file_list <- c("../pool1_nextseq/nextseq2/per_sample_outs/065E-1B-Epsilon-D3/count/sample_filtered_feature_bc_matrix.h5", \
	"../pool1_nextseq/nextseq2/per_sample_outs/065E-1B-Epsilon-D7/count/sample_filtered_feature_bc_matrix.h5", \
	"../pool1_nextseq/nextseq2/per_sample_outs/065E-1B-Epsilon-D10/count/sample_filtered_feature_bc_matrix.h5")

#generate obj individuall using SCT normalization method
seurat_list <- list()
batches <- c("D3", "D7", "D10")
for (i in 1:3) {
  project <- batches[i]
  print(project)
  file <- h5_file_list[i]
  print(file)
  x <- Read10X_h5(file) 
  x <- CreateSeuratObject(x, project = project, min.cells = 3, min.features = 200) %>%
    SCTransform() %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>%
    RunUMAP(dims = 1:30) %>%
    FindClusters()
    RunAzimuth(reference="lungref")

  seurat_list[[project]] <- x
  saveRDS(x, file = paste0(project, "_20230801.rds"))
}

# generate interated obj
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000) 
seurat_list <- PrepSCTIntegration(object.list = seurat_list,
                                  anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = seurat_list,
                                  normalization.method = "SCT",
                                  anchor.features = features)
combined <- IntegrateData(anchorset = anchors)
combined





######## recluster
DefaultAssay(combined) <- "integrated" # switch default assay to the integrated dataset; comments after execution: do not run re-scale using SCT assay
combined <- ScaleData(combined, verbose = FALSE) # if the assay doesn't have scaled data; re-normaize if the assay doesn't have normalized data
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)







#############################################################
# Generate cell counts per types, e.g., cell type, clusters, etc.,
#############################################################
file <- "epsilon_sct_az_hst.rds"
sc <- readRDS(file)
colnames(sc@meta.data)
types <- colnames(sc@meta.data)[c(7,12,14)]

for (i in types){
  print(i)
  cell.counts <- table(sc$orig.ident, sc@meta.data[[i]])
  write.csv(cell.counts, file = paste0("cell_counts_by", i, ".csv"))
}


#############################################################
# using AverageExpression() function to calculate average expression for each cluster
#############################################################
# by default, all the applicable assays including RNA, refAssay, integrated, and etc., will perform the function at the same time
# the returned object e.g., avg.exp is a list of matrix instead of a single matrix
# according to the active idents, and return a matrix for each assay
# if not define "slot = ", the default slot is "data", i.e., log transformed and normalized data
# The ideal assay for calculating AVG gene expression should be the one from "RNA" assay; do not use "integrated" for expression analysis; "SCT" depends.

Idents(d3) <- d3$seurat_clusters
avg.exp <- AverageExpression(d3) # calculaiton AVG on the slot of "data" for all the assays and export a list of matrix
for (assay_id in names(avg.exp)) {
  avg.exp.df <- as.data.frame(avg.exp[[assay_id]])
  write.csv(avg.exp.df, file = paste0("d3", assay_id, "_avg_exp.csv"), row.names = TRUE)
}

avg.exp <- AverageExpression(d3,slot = "counts") # if want to explore using counts slot
for (assay_id in names(avg.exp)) {
  avg.exp.df <- as.data.frame(avg.exp[[assay_id]])
  write.csv(avg.exp.df, file = paste0("d3", assay_id, "_avg_counts.csv"), row.names = TRUE)
}

for (i in batches){
  sc <- seurat_list[[i]]
  Idents(sc) <- sc$predicted.celltype
  avg.exp <- as.data.frame(AverageExpression(sc, assays = "refAssay")) # only peform AVG on 1 assay
  write.csv(avg.exp, file = paste0(run.name, "/", i, "_avg_exp_by_type.csv"))
  }

# This is to perform AVG gene expression on 1 or more conditons
# most important things: define active assay using "RNA" if only 
DefaultyAssay(obj) <- "RNA" # for DE and AVG expression; or define assays = "RNA" inside the function
Idents(sub.obj) <- sub.obj$predicted.celltype
avg.exp <- as.data.frame(AverageExpression(sub.obj, assays = "RNA", group.by = c("days"))) # groups by 1 conditon
write.csv(avg.exp, file = paste0(run.name, "/gene_avg_exp_by_days.csv"))
avg.exp <- as.data.frame(AverageExpression(sub.obj, assays = "RNA", group.by = c("strains")))
write.csv(avg.exp, file = paste0(run.name, "/gene_avg_exp_by_strains.csv"))
avg.exp <- as.data.frame(AverageExpression(sub.obj, assays = "RNA", group.by = c("strains", "days"))) # groups by 2 conditon
write.csv(avg.exp, file = paste0(run.name, "/gene_avg_exp_by_strains_days.csv"))
avg.exp <- as.data.frame(AverageExpression(sub.obj, assays = "RNA", group.by = c("sars.pos")))
write.csv(avg.exp, file = paste0(run.name, "/gene_avg_exp_by_sars.csv"))
avg.exp <- as.data.frame(AverageExpression(sub.obj, assays = "RNA", group.by = c("is.Omicron")))
write.csv(avg.exp, file = paste0(run.name, "/gene_avg_exp_OmicronVSothers.csv"))




##################### CHECK Identical
Omi4MD3 <- seurat_list[["4M-Omicron-D3"]]
Omi4MD3
print(identical(Omi4MD3[["RNA"]]@data, Omi4MD3[["refAssay"]]@data))
print(identical(Omi4MD3[["RNA"]]@counts, Omi4MD3[["refAssay"]]@counts))
print(identical(Omi4MD3[["RNA"]]@data, Omi4MD3[["refAssay"]]@counts))
print(identical(Omi4MD3[["RNA"]]@data, Omi4MD3[["RNA"]]@counts))
print(identical(Omi4MD3[["refAssay"]]@data, Omi4MD3[["refAssay"]]@counts))

[1] FALSE
[1] TRUE
[1] FALSE
[1] FALSE
[1] FALSE

for (i in batches){
  sc <- seurat_list[[i]]
  assay <- sc[["RNA"]]
  print(identical(assay@counts, assay@data))


################# check expression of a list of genes in data
sink(paste0(run.name, "/", "sars_entry_factor.txt"))
print("Expression of the entry factors of Ace2, Temrss2, Furin, Bsg, Nrp1, Ext1 in my dataset")
for (i in batches){
  sc <- seurat_list[[i]]
  print(i)
  numbers <- features %in% rownames(sc[["RNA"]])
  print(numbers)
  print(sum(numbers))
  print("")
}
sink()

sars <- c("S", "M", "N", "ORF1ab", "ORF8")
sink("sars_check.txt")
print("check the expression of sars genes S, M, N, ORF1ab, and ORF8 in my datasets")
for (i in batches){
  sc <- seurat_list[[i]]
  print(i)
  print("dimension of refAssay")
  print(dim(sc[["refAssay"]]))
  print("dimension of RNA")
  print(dim(sc[["RNA"]]))
  numbers <- sars %in% rownames(sc[["RNA"]])
  print(numbers)
  print(sum(numbers))
  print("")
}
sink()

# save the output of a block to a file
# note the behaviors of "sink": it only save the stdout of function execution, but not save the desired outputs files/figs/folder from functions
# do not use "sink" if the functions will export files
sink("folder/filename.txt")
print("a comment")
table(obj@meta.data$strains)
sink()



## downsampleing
1.random sampleing only define number
downsampled.obj <- large.obj[, sample(colnames(large.obj), size = ncol(small.obj), replace=F)]

# downsampleing to make sure each categroy has the same number
# subset meta dataframe; each row of the meta is a cell, each column is its attribute, e.g., grouping, celltype, clusteringinfo, etc.,
# subset meta dataframe by condition  for D10, get its cell id
d10 <- rownames(meta[meta$days == "D10",])
# subset d3 cells's meta, do sampleing based on total number of d10 cells, and get its cell id
d3 <- meta[meta$days == "D03",] %>% 
  sample_n(size = length(d10), replace = FALSE) %>%
  rownames()
print(length(d3))

# subset d7 cells's meta, do sampleing based on total number of d10 cells, and get its cell id
d7 <- meta[meta$days == "D07",] %>% 
  sample_n(size = length(d10), replace = FALSE) %>%
  rownames()
print(length(d7))

# subset the whole seurat object using the above retrieved cell id
to.keep <- c(d3, d7, d10) # cell id to keep
sub.obj2 <- subset(sub.obj, cells = to.keep)


# using leiden
sinteractive --mem=100g --cpus-per-task=12 --gres=lscratch 100
cd $IRF
source myconda
conda active r4seurat
module load rstudio R
rstudio&

library(reticulate)
use_condaenv("/data/pengl7/conda/envs/r4seurat")
py_config()
# need to first create a conda env as r4seurat; conda active r4seurat; conda install leidenalg; conda install numpy pandas; start rstudio&
leidenalg <- import("leidenalg")


# remove unless objs
rm(virus.seurat.list)

# integrate extra large datasets (refer to 20230831-ingetrate: integrate_virsu.rmd)
https://satijalab.org/seurat/articles/integration_large_datasets.html
Since this dataset contains D3, D7, D10, we will chose one from each groups to use in a reference-based workflow; 
this is a faster otpiton for very large datset; introduce here the possibility of specifying one or more of the datasets as the ‘reference’ for integrated analysis, 
with the remainder designated as ‘query’ datasets. 
For example, when integrating 10 datasets with one specified as a reference, we perform only 9 comparisons. 
my case 10 samples of control are not that big; if I still chose to std wf which identify anchors between all pairs of datasets; 
we will perform 45 different pairwise comparisons

# need to run SelectIntegrationFeatures and FindIntegrationAnchors
virus.seurat.list <- list()
for (i in virus){
  file <- i
  id <- gsub(pattern = "^2023-08-28-|^2023-08-29-|\\.rds$", "", i)
  x <- readRDS(paste0("../20230823-individual/individuals-2/", i))
  DefaultAssay(x) <- "SCT" # have to switch to SCT so that no error in IntegrateData()
  virus.seurat.list[[id]] <- x
}
features <- SelectIntegrationFeatures(object.list = virus.seurat.list)
virus.seurat.list <- lapply(X = virus.seurat.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = virus.seurat.list, reference = c(1, 2, 3), reduction = "rpca", dims = 1:50)
obj.integrated <- IntegrateData(anchorset = anchors, dims = 1:50, normalization.method = "SCT") # taks very long time for IntegrateData
obj.integrated <- ScaleData(obj.integrated, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:50)

# Perform DE analysis using SCT # not sure whether I used it successfuly or not
https://satijalab.org/seurat/reference/prepsctfindmarkers
NEED PrepSCTFindMarkers()
Idents(sub.obj) <- sub.obj$days
sub.obj <- PrepSCTFindMarkers(object = sub.obj) # takes about 1 h


# plot histogram
hist(sars.avg)

# Add the new column to the meta.data DataFrame
sub.obj@meta.data[["sars.avg"]] <- sars.avg
sub.obj@meta.data[["sars.pos"]] <- ifelse(sub.obj$sars.avg > 0, "POS", "NEG")
sub.obj@meta.data[["is.Omicron"]] <- ifelse(sub.obj$strains != "Omicron", "Others", sub.obj$strains)


# some math
sars.avg <- colMeans(sub.obj[["RNA"]][sars,])
length(sars.avg)
range(sars.avg)
quantile(sars.avg, probs = 0.60)

# check filesize
format(object.size(integrated), units = "Gb")