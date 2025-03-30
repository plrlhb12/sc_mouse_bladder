# arrange multiple plots in Seurat using patchwork
library(patchwork)

p1 <- FeaturePlot(pbmcsca, features = "CCR7")
p2 <- FeaturePlot(pbmcsca, features = "FCGR3A")
p3 <- VlnPlot(pbmcsca, features = "AXL", group.by = "predicted.celltype.l2", idents = c("ASDC",
    "pDC", "cDC1", "cDC2"))
p4 <- FeaturePlot(pbmcsca, features = "predictionscorecelltypel2_Treg")

p1 + p2 + p3 + p4 + plot_layout(ncol = 2)
(p1 + p2) & NoLegend()
(p1 + p2)/p3