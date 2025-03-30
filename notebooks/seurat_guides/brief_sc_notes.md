# Our single cell gene expression application is 10x Genomics probe-based fixed RNA profiling, so mitochondrial genes doesn't exist in the profiled gene list. No need to filter cells based on the expression percentage of mitochondrial genes
# I added the SARS-CoV-2 genes (M, N, S, ORF1ab, ORF8) into the feature table.

# 1. generate and save an individual seurat object for each sample using default or minimal QC thresholds

`x <- Read10X_h5(file_path) 
  x <- CreateSeuratObject(x, project = project, min.cells = 3, min.features = 200) %>%
    SCTransform() %>%
    RunPCA() %>%
    FindNeighbors(dims = 1:30) %>%
    RunUMAP(dims = 1:30) %>%
    FindClusters(algorithm = 4, selection.method = "vst", nfeatures = 2000) %>%  # using leiden clustering
    RunAzimuth(reference="../hamsterLungRef") # perform annotation using a reference object from a paper
`
# 2. read and integrate all SARS rds as the "virus" obj, add meta information (e.g., strains and expression status (POS or NEG) of SARS genes), and performed differential gene expression analysis

`virus.list <- list()
  for (i in virus){
    file <- i
    id <- gsub(pattern = "^2023-08-28-|^2023-08-29-|\\.rds$", "", i)
    x <- readRDS(paste0("../20230823-individual/individuals-2/", i))
    DefaultAssay(x) <- "SCT"
    virus.seurat.list[[id]] <- x
  }
  anchors <- FindIntegrationAnchors(object.list = virus.list, reference = c(1, 2, 3), reduction = "rpca", dims = 1:50)
  virus.integrated <- IntegrateData(anchorset = anchors, dims = 1:50, normalization.method = "SCT")
`
# 3. read and integrate all uninfected rds as the "ctrl" obj and performed differential gene expression analysis

`anchors <- FindIntegrationAnchors(object.list = ctrl.list, reference = c(1, 2, 3), reduction = "rpca", dims = 1:50)
  ctrl.integrated <- IntegrateData(anchorset = anchors, dims = 1:50, normalization.method = "SCT")
`
# 4. integrate the above ctrl and virus rds objects as the "all_cleaned" and peform differential gene expression

  `virus.seurat.list <- c(ctrl, virus)
  features <- SelectIntegrationFeatures(object.list = virus.seurat.list)
  anchors <- FindIntegrationAnchors(object.list = virus.seurat.list, reduction = "rpca", dims = 1:50)
  virus.seurat.list <- lapply(X = virus.seurat.list, FUN = function(x) {
    DefaultAssay(x) <- "SCT"
  })
  all_cleaned <- IntegrateData(anchorset = anchors, dims = 1:50)
  `


