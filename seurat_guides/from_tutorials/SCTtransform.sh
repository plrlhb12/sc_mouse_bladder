# Performing integration on datasets normalized with SCTransform
# SCTransfrom is an improved method for the normalization of scRNA-seq based on regularized negative binomial regression
# It avoids some of the pitfalls of standard normalization workflows, including the addition of a pseudocount, and log-transformation
# This procedure omits the need for heuristic steps including pseudocount addition or log-transformation 
# and improves common downstream analytical tasks such as variable gene selection, dimensional reduction, and differential expression.
# https://satijalab.org/seurat/articles/sctransform_vignette.html


############ general example using v1

library(Seurat)
library(ggplot2)
library(sctransform)

pbmc_data <- Read10X(data.dir = "../data/pbmc3k/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc_data)
# store mitochondrial percentage in object meta data
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")
# run sctransform
# Note that this single command replaces NormalizeData(), ScaleData(), and FindVariableFeatures().
# Transformed data will be available in the SCT assay; also can remove confounding sources of variation, for example, mitochondrial mapping percentage
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
# using glmGamPoi package which substantially improves the speed of the learning procedure
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("glmGamPoi")
pbmc <- SCTransform(pbmc, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
# These are now standard steps in the Seurat workflow for visualization and clustering
pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE)
DimPlot(pbmc, label = TRUE) + NoLegend()

# in a single command
pbmc <- CreateSeuratObject(pbmc_data) %>%
    PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
    SCTransform(vars.to.regress = "percent.mt") %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>%
    RunUMAP(dims = 1:30) %>%
    FindClusters()

############ example using v2 regularization which improves compared v1
# install sctransform >= 0.3.3
install.packages("sctransform")
# invoke sctransform - requires Seurat>=4.1
object <- SCTransform(object, vst.flavor = "v2")
# also install the glmGamPoi package which substantially improves the speed of the learning procedure.
$$$$$$$$$$$$$$$$$$$$$$$ not copying the whole tutorial yet
https://satijalab.org/seurat/articles/sctransform_v2_vignette.html

############ example in integration: modify Seurat integration workflow with SCTransfromed function

LoadData("ifnb")
ifnb.list <- SplitObject(ifnb, split.by = "stim")
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
    anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30)
p1 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "seurat_annotations", label = TRUE,
    repel = TRUE)
p1 + p2


###########