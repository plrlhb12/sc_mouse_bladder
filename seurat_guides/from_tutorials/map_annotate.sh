# Mapping and annotating query datasets
# need high quality reference dataset to integrate and be mapped to
# doesn't need correct underlying raw query data
# Steps: 
# 1. generate an integrated reference
# 2. cell type label transfer


# example: human pancratic islet datsets across 4 technologies

################### Part 1: data integration
################### if the reference dataset is a single dataset, no need to do split and integration later
# This example split 1 dataset into 4; of which 3 datasets are integrate into 1 reference data

library(Seurat)
library(SeuratData)
InstallData("panc8")
data("panc8")

# generate a list of datasets for integration
pancreas.list <- SplitObject(panc8, split.by = "tech")
pancreas.list <- pancreas.list[c("celseq", "celseq2", "fluidigmc1", "smartseq2")]
reference.list <- pancreas.list[c("celseq", "celseq2", "smartseq2")] # need to integrate as 1 reference
pancreas.query <- pancreas.list[["fluidigmc1"]] # this is the query dataset

# perfrom std normalization and identify HVG
for (i in 1:length(pancreas.list)) {
    pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
    pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", nfeatures = 2000,
        verbose = FALSE)
}

############ generate the reference dataset having celltype info, for future real single reference dataset, no need integration
# first integrate 3 dataset as 1 reference dataset, apply std sc pipeline
# notice that this sample dataset already contain "celltype" metadata

# for any integration in Seurat, need to use FindIntegrationAnchor() and IntegrateData()
# features <- SelectIntegrationFeatures(object.list = ifnb.list) this was not included in this tutorial
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30) # default dim is 30, can range from 10 to 50
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)
library(ggplot2)
library(cowplot)
library(patchwork)
DefaultAssay(pancreas.integrated) <- "integrated"
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) +
    NoLegend()
p1 + p2


####################### Data transfer or label transfere
# notice the distinctions between data transfer and data integration:
# 1. data integration: FindTransferAnchor() and IntegrateData(), correct dataset in the 2nd funciton
# 2. data/label transfer: FindTransferAnchor() and TransferData(), then AddMetaData(); Seurat does not correct or modify the query expression data.
	# Seurat default to project the PCA structure of a reference onto the query, instead of learning a joint structure with CCA. We generally suggest using this option when projecting data between scRNA-seq datasets.

pancreas.anchors <- FindTransferAnchors(reference = pancreas.integrated, query = pancreas.query,
    dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = pancreas.anchors, refdata = pancreas.integrated$celltype,
    dims = 1:30)
pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)
pancreas.query$prediction.match <- pancreas.query$predicted.id == pancreas.query$celltype
table(pancreas.query$prediction.match)
# verify further by visual canonical cell type markers
#  Note that even though some of these cell types are only represented by one or two cells (e.g. epsilon cells), we are still able to classify them correctly
table(pancreas.query$predicted.id)
VlnPlot(pancreas.query, c("REG1A", "PPY", "SST", "GHRL", "VWF", "SOX10"), group.by = "predicted.id")

##################### UMAP: Unimodal UMAP Projection: project query onto the reference UMAP structure
# generate UMAP on reference
pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:30, reduction = "pca", return.model = TRUE)
# map query to reference UMAP
pancreas.query <- MapQuery(anchorset = pancreas.anchors, reference = pancreas.integrated, query = pancreas.query,
    refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap")
# now can visualize query cells alongside reference 
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 3,
    repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(pancreas.query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
    label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2