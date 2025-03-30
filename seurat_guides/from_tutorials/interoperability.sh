# install scater https://bioconductor.org/packages/release/bioc/html/scater.html
library(scater)
library(Seurat)
# install SeuratDisk from GitHub using the remotes package remotes::install_github(repo = 'mojaveazure/seurat-disk', ref = 'develop')
library(SeuratDisk)
library(SeuratData)
library(patchwork)

################# Seurat vs SingleCellExperiment
# converting to/from SingleCellExperiment using as.SingleCellExperiment()
# Use PBMC3K from SeuratData and conver it to SingelCellExperiment for use with scater package
# which focus on with a focus on quality control and visualization of sc RNA-seq data
InstallData("pbmc3k")
pbmc <- LoadData(ds = "pbmc3k", type = "pbmc3k.final")
pbmc.sce <- as.SingleCellExperiment(pbmc)
p1 <- plotExpression(pbmc.sce, features = "MS4A1", x = "ident") + theme(axis.text.x = element_text(angle = 45,
    hjust = 1))
p2 <- plotPCA(pbmc.sce, colour_by = "ident")
p1 + p2


# convert SingleCellExperiment to Seurat using as.Seurat()
# download from hemberg lab https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/manno_human.rds
manno <- readRDS(file = "../data/manno_human.rds")
manno <- runPCA(manno)
manno.seurat <- as.Seurat(manno, counts = "counts", data = "logcounts")
# gives the same results; but omits defaults provided in the last line
manno.seurat <- as.Seurat(manno)
Idents(manno.seurat) <- "cell_type1"
p1 <- DimPlot(manno.seurat, reduction = "PCA", group.by = "Source") + NoLegend()
p2 <- RidgePlot(manno.seurat, features = "ACTB", group.by = "Source")
p1 + p2


##################  Converting to/from loom using as.loom
# to loom
pbmc.loom <- as.loom(pbmc, filename = "../output/pbmc3k.loom", verbose = FALSE)
pbmc.loom
pbmc.loom$close_all() # Always remember to close loom files when done

# from loom: download from linnarsson lab using as.Seurat
# https://storage.googleapis.com/linnarsson-lab-loom/l6_r1_immune_cells.loom
l6.immune <- Connect(filename = "../data/l6_r1_immune_cells.loom", mode = "r")
l6.immune
l6.seurat <- as.Seurat(l6.immune)
Idents(l6.seurat) <- "ClusterName"
VlnPlot(l6.seurat, features = c("Sparc", "Ftl1", "Junb", "Ccl4"), ncol = 2)
l6.immune$close_all() # Always remember to close loom files when done

#################  Converting to/from AnnData
# https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html
library(Seurat)
library(SeuratData)
library(SeuratDisk)

# From Seurat to Anndata, need 2 steps. 1st save Seurat as an h5Seurat, 2nd convert it to AnnData
SaveH5Seurat(pbmc3k.final, filename = "pbmc3k.h5Seurat") # save a file
Convert("pbmc3k.h5Seurat", dest = "h5ad")

import scanpy # in python
adata = scanpy.read_h5ad("pbmc3k.h5Seurat")

# From Anndata to Seurat, also 2 steps. First convert Anndata h5ad to h5Seurat then LoadH5Seurat
Convert("pbmc3k_final.h5ad", dest = "h5seurat", overwrite = TRUE)
pbmc3k <- LoadH5Seurat("pbmc3k_final.h5seurat")

