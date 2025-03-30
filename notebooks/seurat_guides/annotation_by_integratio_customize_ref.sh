#####################################################
# annotation method 1 using Azimuth
# https://satijalab.github.io/azimuth/articles/run_azimuth_tutorial.html
# using built-in ref
lung72 <- RunAzimuth(lung72, reference = "lungref")

p1 <- DimPlot(lung72, group.by = "predicted.ann_level_1", label = TRUE, label.size = 3) + NoLegend()
p2 <- DimPlot(lung72, group.by = "predicted.ann_level_2", label = TRUE, label.size = 3) + NoLegend()
p3 <- DimPlot(lung72, group.by = "predicted.ann_level_3", label = TRUE, label.size = 3) + NoLegend()
p4 <- DimPlot(lung72, group.by = "predicted.ann_level_4", label = TRUE, label.size = 3) + NoLegend()
p1 + p2
p3 + p4
p5 <- DimPlot(lung72, group.by = "predicted.ann_level_5", label = TRUE, label.size = 3) + NoLegend()
p5

saveRDS(object = lung72, file = "lung72h-AZ")


# using customize ref
check my samples of hamsterRef


####################################################
# annotation method 2: using integration
# https://satijalab.org/seurat/articles/integration_mapping.html
library(lungref.SeuratData)
lung_ref <- LoadData("lungref", "azimuth")
lung_ref <- NormalizeData(lung_ref$map)
colnames(lung_ref[["map"]]@meta.data)
lung.anchors <- FindTransferAnchors(reference = lung_ref$map, query = lung72,
                                    dims = 1:30, reference.reduction = "pca")
