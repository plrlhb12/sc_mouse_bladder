################## ################## ################## ##################
# This function work with v4
# Azimuth: has reference data of RNA,ADT, ATAC bride for human PBMC, Motor Cortex, Pancreas, lung, kidney, liver, heart, tonsil, bone marrow and fetal
# using Reference-based mapping 

################### two ways, app and cmd; app has limitation on cell numbers (10,000)
# https://azimuth.hubmapconsortium.org/references/                         # app version
# https://satijalab.github.io/azimuth/articles/run_azimuth_tutorial.html    # cmd version
# 
# for users interested in performing these analyses outside the context of the Azimuth app, we suggest using Seurat v4 and using our vignette on Mapping and annotating query datasets as an example. 
# You can also download a Seurat v4 R script on the “Download Results” from the app once your analysis is complete to reproduce the results locally.

# We do not recommend pre-filtering genes in the data you upload. Users should upload an unprocessed counts matrix, 
# for example, the output of the 10x Genomics CellRanger pipeline
# must have unnormalized "count" data in the "RNA" assay slot

# Cell types must have at least 15 query cells of that predicted type to find biomarkers.

################## cmd version, using only one line cmd: obj <- RunAzimuth(obj, reference = "")
# for sc RNA-Seq, reference datasets are automatically download from SeuratData framework; work with V5 and V5
# for atac, transfer annotation from RNA references to ATAC query; only work with V5

devtools::install_github("satijalab/seurat", "seurat5")
devtools::install_github("satijalab/seurat-data", "seurat5")
devtools::install_github("satijalab/azimuth", "seurat5")
library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)

# generate or load an Seurat object, even direct take the path of h5ad as input
# keep in mind that need unnormalized counts from the object
bm <- RunAzimuth(query = "~/human_cd34_bone_marrow.h5ad", reference = "bonemarrowref")
pbmcsca <- RunAzimuth(pbmcsca, reference = "pbmcref")

# check available references
##                                  Dataset Version                        Summary
## adiposeref.SeuratData         adiposeref   1.0.0     Azimuth Reference: adipose
## bonemarrowref.SeuratData   bonemarrowref   1.0.0  Azimuth Reference: bonemarrow
## fetusref.SeuratData             fetusref   1.0.0       Azimuth Reference: fetus
## heartref.SeuratData             heartref   1.0.0       Azimuth Reference: heart
## humancortexref.SeuratData humancortexref   1.0.0 Azimuth Reference: humancortex
## kidneyref.SeuratData           kidneyref   1.0.1      Azimuth Reference: kidney
## lungref.SeuratData               lungref   2.0.0        Azimuth Reference: lung
## mousecortexref.SeuratData mousecortexref   1.0.0 Azimuth Reference: mousecortex
## pancreasref.SeuratData       pancreasref   1.0.0    Azimuth Reference: pancreas
## pbmcref.SeuratData               pbmcref   1.0.0        Azimuth Reference: pbmc
## tonsilref.SeuratData           tonsilref   1.0.0      Azimuth Reference: tonsil

available_data <- AvailableData()
available_data[grep("Azimuth", available_data[, 3]), 1:3]

# check the numbers of each predition cluster
sort(table(bm$predicted.celltype.l2), decreasing = TRUE)
# check predition score
bm$predicted.celltype.l2.scores

# visualize; remember to normalize before any visualization
pbmcsca <- RunAzimuth(pbmcsca, reference = "pbmcref")
pbmcsca <- NormalizeData(pbmcsca)
Idents(pbmcsca) <- "predicted.celltype.l2"   # change Idents
p1 <- FeaturePlot(pbmcsca, features = "CCR7")
p2 <- FeaturePlot(pbmcsca, features = "FCGR3A")
p3 <- VlnPlot(pbmcsca, features = "AXL", group.by = "predicted.celltype.l2", idents = c("ASDC",
    "pDC", "cDC1", "cDC2"))
p4 <- FeaturePlot(pbmcsca, features = "predictionscorecelltypel2_Treg")
p1 + p2 + p3 + p4 + plot_layout(ncol = 2)

p1 <- DimPlot(pbmcsca, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()
p2 <- DimPlot(pbmcsca, group.by = "Method")
p1 + p2

################### FAQ
# Should I map my batches separately or combined?
We have observed that the results (both for visualization and annotation) are very similar when mapping individual batches separately, or combined together
Azimuth can successfully remove batch effects between query and reference cells, even when there are multiple query batches. However, as discussed further below, the results of QC metrics may change

# What optimizations are in the app that are not default in Seurat?
To optimize the web app time and resource consumption, we made several changes to the base Seurat mapping workflow.

Decrease genes and cells used for SCTtranform
Downsample reference and decrese n.trees; downsample query dataset to at most 5,000 cells in app 
Use presto package for DE, which analyze at least 15 query cells


# What if my query dataset contains cell types that aren't present in the reference?
the mapping will likely return poor ‘Dataset-level’ QC metrics including low mapping scores

# How can I tell if my mapping results are accurate?
a single metric is insufficient to describe the quality of mapping
users should not limit their evaluation of mapping to these QC metrics, and can should explore their results

Dataset-level Metrics: 
1. % of query cells with anchors:Typically, we observe values >15% when mapping is successful.  
However, in some cases, when there is a large batch effect between query and reference datasets from the same tissue, this metric can fall below 15% even when mapping is successful.
2. Cluster preservation score: Scores are scaled from 0 (poor) to 5 (best).
If the query dataset consists of a homogeneous group of cells, or if the query dataset contains cells from multiple batches (which would be corrected by Azimuth), this metric may return a low value even in cases where mapping is successful. 

The score is calculated using the ClusterPreservationScore function in Azimuth.
Cell-level Metrics:
1. Prediction scores: Cell prediction scores range from 0 to 1 and reflect the confidence associated with each annotation. Cells with high-confidence annotations (for example, prediction scores > 0.75) reflect predictions that are supported by mulitple consistent anchors.
2. Mapping scores: This value from 0 to 1 reflects confidence that this cell is well represented by the reference.

How can a cell get a low prediction score and a high mapping score?
A cell can get a low prediction score because its probability is equally split between two clusters (for example, for some cells, it may not be possible to confidently classify them between the two possibilities of CD4 Central Memory (CM), and Effector Memory (EM), which lowers the prediction score, but the mapping score will remain high.

How can a cell get a high prediction score and a low mapping score?
A high prediction score means that a high proportion of reference cells near a query cell have the same label. However, these reference cells may not represent the query cell well, resulting in a low mapping score. Cell types that are not present in the reference should have lower mapping scores. For example, we have observed that query datasets containing neutrophils (which are not present in our reference), will be confidently annotated as CD14 Monocytes, as Monocytes are the closest cell type to neutrophils, but receive a low mapping score.





################## ################## ################## ##################



