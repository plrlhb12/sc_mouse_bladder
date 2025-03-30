# Data visualization methods in Seurat

library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
SeuratData::InstallData("pbmc3k")
data("pbmc3k.final")
pbmc3k.final$groups <- sample(c("group1", "group2"), size = ncol(pbmc3k.final), replace = TRUE)
features <- c("LYZ", "CCL5", "IL32", "PTPRCAP", "FCGR3A", "PF4")
pbmc3k.final
############################# Tips in chunks of Rstudio
1. ggsave comes from pakage of ggplot2
2. the dim of the plot shown in console is different from that saved in local, so often need to ajust width and height when save the plots
3. the plot in console may not detailed enough or distoreted, can zoom it out to see the whole picture
4. ggsave take the current dev in the current chunk as input, so 

FeaturePlot(combined, features = sars, split.by = "orig.ident")
ggsave(filename = "sars_by_days.pdf", width = 6, height = 12)
#is the same as 
p1 <- FeaturePlot(combined, features = sars, split.by = "orig.ident")
ggsave(filename = "sars_by_days.pdf", p1, width = 6, height = 12)

# or
p1 <- DimPlot(combined, reduction = "umap", group.by = "orig.ident") # umap, color by condtion
p2 <- DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE) # umap, color by clusters
p1 + p2
ggsave("episilon_umap.pdf", width = 10, height = 4)

# arrange multiple plots
p1 <- DimPlot(combined, reduction = "umap", group.by = "orig.ident") # umap, color by condtion
p2 <- DimPlot(combined, reduction = "umap", repel = TRUE, group.by = "seurat_clusters") # umap, color by clusters
p3 <- DimPlot(combined, reduction = "umap", group.by = "predicted.ann_level_2")
p4 <- DimPlot(combined, reduction = "umap", group.by = "predicted.ann_level_3")
a_plot <- p1 + p2 + p3 + p4 + plot_layout(ncol = 2)
ggsave("episilon_umap_az.pdf", a_plot, width = 16, height = 12)


p1 <- DimPlot(d3, label = TRUE)
p2 <- DimPlot(d3, label = TRUE, group.by = "predicted.celltype", label.size = 2)
p1 + p2
ggsave("4J-Omicron-D3/20230818/18cluster_umap.pdf",width = 15) + theme(text = element_text(size = 8))

#############################  Five visualizations of marker feature expression
# Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
RidgePlot(pbmc3k.final, features = features, ncol = 2)

# Violin plot - Visualize single cell expression distributions in each cluster
VlnPlot(pbmc3k.final, features = features)
# Violin plots can also be split on some variable. Simply add the splitting variable to object metadata and pass it to the split.by argument
VlnPlot(pbmc3k.final, features = "percent.mt", split.by = "groups")

####### FeaturePlot
# Feature plot - visualize feature expression in low-dimensional space, e.g., UMAP
FeaturePlot(pbmc3k.final, features = features)
FeaturePlot(pbmc3k.final, features = "MS4A1", min.cutoff = 1, max.cutoff = 3) # Adjust the contrast in the plot
# Calculate feature-specific contrast levels based on quantiles of non-zero expression. Particularly useful when plotting multiple markers
FeaturePlot(pbmc3k.final, features = c("MS4A1", "PTPRCAP"), min.cutoff = "q10", max.cutoff = "q90")
# Visualize co-expression of two features simultaneously
FeaturePlot(pbmc3k.final, features = c("MS4A1", "CD79A"), blend = TRUE)
# Split visualization to view expression by groups (replaces FeatureHeatmap)
FeaturePlot(pbmc3k.final, features = c("MS4A1", "CD79A"), split.by = "groups")

# Dot plots - the size of the dot corresponds to the percentage of cells expressing the feature in each cluster. The color represents the average expression level
DotPlot(pbmc3k.final, features = features) + RotatedAxis()
# SplitDotPlotGG has been replaced with the `split.by` parameter for DotPlot
DotPlot(pbmc3k.final, features = features, split.by = "groups") + RotatedAxis()

# Single cell heatmap of feature expression
DoHeatmap(subset(pbmc3k.final, downsample = 100), features = features, size = 3)
# DoHeatmap now shows a grouping bar, splitting the heatmap into groups or clusters. This can be changed with the `group.by` parameter
DoHeatmap(pbmc3k.final, features = VariableFeatures(pbmc3k.final)[1:100], cells = 1:500, size = 4, angle = 90) + NoLegend()


# DimPlot replaces TSNEPlot, PCAPlot, etc. In addition, it will plot either 'umap', 'tsne', or 'pca' by default, in that order
DimPlot(pbmc3k.final) # default plot umap

pbmc3k.final.no.umap <- pbmc3k.final
pbmc3k.final.no.umap[["umap"]] <- NULL # if no umap, no tsne, plot pca
DimPlot(pbmc3k.final.no.umap) + RotatedAxis()


############################ plotting multiple plots together using patchwork system instead of previous CombinePlot()
plot1 <- DimPlot(pbmc3k.final)
plot2 <- FeatureScatter(pbmc3k.final, feature1 = "LYZ", feature2 = "CCL5")
plot1 + plot2 # Combine two plots
(plot1 + plot2) & NoLegend() # Remove the legend from all plots

############################ plotting accessories
# LabelClusters and LabelPoints will label clusters (a coloring variable) or individual points on a ggplot2-based scatter plot
plot <- DimPlot(pbmc3k.final, reduction = "pca") + NoLegend()
LabelClusters(plot = plot, id = "ident")

# Both functions support `repel`, which will intelligently stagger labels and draw connecting lines from the labels to the points or clusters
LabelPoints(plot = plot, points = TopCells(object = pbmc3k.final[["pca"]]), repel = TRUE)


############################ adding themes to plots
baseplot <- DimPlot(pbmc3k.final, reduction = "umap")
# Add custom labels and titles
baseplot + labs(title = "Clustering of 2,700 PBMCs")
# Use community-created themes, overwriting the default Seurat-applied theme Install ggmin with remotes::install_github('sjessa/ggmin')
baseplot + ggmin::theme_powerpoint()
# Seurat also provides several built-in themes, such as DarkTheme; for more details see ?SeuratTheme
baseplot + DarkTheme()
# Chain themes together
baseplot + FontSize(x.title = 20, y.title = 20) + NoLegend()


################  interactive plotting
# Include additional data to display alongside cell names by passing in a data frame of information Works well when using FetchData
plot <- FeaturePlot(pbmc3k.final, features = "MS4A1")
HoverLocator(plot = plot, information = FetchData(pbmc3k.final, vars = c("ident", "PC_1", "nFeature_RNA")))

pbmc3k.final <- RenameIdents(pbmc3k.final, DC = "CD14+ Mono")
plot <- DimPlot(pbmc3k.final, reduction = "umap")
select.cells <- CellSelector(plot = plot)
................ 


################
p3 <- VlnPlot(pbmcsca, features = "AXL", group.by = "predicted.celltype.l2", idents = c("ASDC",
    "pDC", "cDC1", "cDC2"))
p4 <- FeaturePlot(pbmcsca, features = "predictionscorecelltypel2_Treg")
p1 + p2 + p3 + p4 + plot_layout(ncol = 2)

p1 <- DimPlot(pbmcsca, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()
p2 <- DimPlot(pbmcsca, group.by = "Method")
p1 + p2


############## visualize a specific cell type or cluster
# since the color are so simlliar, show the types separately
for (type in levels(factor(d3$predicted.celltype))){
  DimPlot(d3, group.by = "predicted.celltype", cells = d3$predicted.celltype == type) +xlim(-15,15) +ylim(-15,15)
  ggsave(paste0("4J-Omicron-D3/20230818/", type, "_umap.pdf"))
}


# DoHeatmap need scaled data
#d3 <- ScaleData(d3)
DefaultAssay(d3) <- "SCT" # SCT can be used for visualization, however be cautions if there are seveal dataset together with big difference of sequencing depth
DoHeatmap(d3, features = top10_cluster_markers$gene, size = 5)+NoLegend() + theme(text = element_text(size = 8))
ggsave(filename = "4J-Omicron-D3/20230818/4J-Omicron-D3-heatmap_top10_18clusters.pdf", height = 14, width = 12)


# batch process
for (i in batches){
  sc <- seurat_list[[i]]
  DefaultAssay(sc) <- "RNA"
  p1 <- VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA"), group.by = "orig.ident", ncol = 2)
  ggsave(filename = paste0(run.name, "/", i, "_QC_violin.pdf"), plot = p1, width = 4, height = 4)
  p2 <- FeatureScatter(sc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  ggsave(filename = paste0(run.name, "/", i, "_QC_scatter.pdf"), plot = p2, width = 4, height = 4)
  p4 <- ElbowPlot(sc)
  ggsave(filename = paste0(run.name, "/", i,"_elbow.pdf"), plot = p4, width = 4, height = 4)
  p5 <- DimPlot(sc, label = TRUE)
  p6 <- DimPlot(sc, label = TRUE, group.by = "predicted.celltype", label.size = 2)
  ggsave(filename = paste0(run.name, "/", i,"_umap.pdf"), plot = (p5 + p6), width = 15) + theme(text = element_text(size = 8))
}


features <- c("S", "M", "N", "ORF1ab", "ORF8")
for (i in batches){
  sc <- seurat_list[[i]]
  DefaultAssay(sc) <- "RNA"
  p1 <- RidgePlot(sc, features = features)
  ggsave(filename = paste0(run.name, "/", i,"_sars_ridge.pdf"), plot = p3, height = 8)
  p2 <- FeaturePlot(sc, features = features)
  ggsave(filename = paste0(run.name, "/", i,"_sars.pdf"), plot = p7)
  p3 <- VlnPlot(sc, features = features)
  ggsave(filename = paste0(run.name, "/", i,"_sars_violin.pdf"), plot = p3, height = 8)
}

