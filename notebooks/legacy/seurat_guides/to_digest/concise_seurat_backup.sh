# see the functions of Seurat
# https://satijalab.org/seurat/reference/index.html


# input
Read10X() # for output dir from Cellranger
Read10X_h5() # for .h5 fle from Cellranger
readRDS() # for .rds saved an R object from Seurat
pbmc3k <- LoadData("pbmc3k") # from build-in dataset
hamsterlung <- readRDS("../../seurat_biowulf/FixedRNAPilot_seurat/lung72h.rds")
RunAzimuth(file.h5ad, reference="lungref")

# std
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))




# incorporate SCTransform









# integration






# label transfer






