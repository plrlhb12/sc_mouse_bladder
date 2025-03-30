
#1.random sampleing only define number
downsampled.obj <- large.obj[, sample(colnames(large.obj), size = ncol(small.obj), replace=F)]

#2. real cases to make sure each category has the same number of cells

#2.1 downsample days to make sure each day has the same number of cells
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


#2.2 downsample sars to make sure POS and NEG have the same mount of cells

# downsample sars POS
pos <- rownames(meta[meta$sars.pos == "POS", ]) 
neg <- meta[meta$sars.pos == "NEG", ] %>%
  sample_n(size = length(pos), replace = FALSE) %>%
  rownames()
to.keep <- c(neg, pos)

down.obj <- subset(sub.obj, cells = to.keep)
DimPlot(down.obj, pt.size = 1, group.by = "predicted.celltype", split.by = "sars.pos")
ggsave(filename = "all_virus/downsample_sars_umap.pdf", width = 9)


#2.2 downsample sars to make sure POS and NEG have the same mount of cells
# downsample for Omicron to make sure Omicron and non-Omicron have the same number of cells
pos <- rownames(meta[meta$is.Omicron == "Omicron", ]) 
neg <- meta[meta$is.Omicron == "Others", ] %>%
  sample_n(size = length(pos), replace = FALSE) %>%
  rownames()
to.keep <- c(neg, pos)

down.obj <- subset(sub.obj, cells = to.keep)
DimPlot(down.obj, group.by = "predicted.celltype", split.by = "is.Omicron")
ggsave(filename = "all_virus/downsample_OmicronVsOthers_umap.pdf", width = 9)




