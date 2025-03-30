#  **Most important things on executing Seurat**
## 1. **using RNA assay for gene expression and its comparison AND visualization**
- for DE
	DefaultyAssay(obj) <- "RNA"

- for AVG gene expression (a little bit more trick): need to specific define the assay inside the function
	DefaultyAssay(obj) <- "RNA" # this doesn't help because the function of AverageExpression() will perform calculation on all available assays in a seurat object

	option 1: define the assay inside the function
	```
	# default is to perform analysis on all assays on "data" slot
	avg.exp <- as.data.frame(AverageExpression(sub.obj, assays = "RNA", group.by = c("strains", "days"))) # groups by 2 conditon
	write.csv(avg.exp, file = paste0(run.name, "/gene_avg_exp_by_strains_days.csv"))
	```

	option 2: subset the output of AVG expression
	```
	Idents(d3) <- d3$seurat_clusters
	avg.exp <- AverageExpression(d3) # calculaiton AVG on the slot of "data" for all the assays and export a list of matrix
	for (assay_id in names(avg.exp)) {
  	avg.exp.df <- as.data.frame(avg.exp[[assay_id]])
  	write.csv(avg.exp.df, file = paste0("d3", assay_id, "_avg_exp.csv"), row.names = TRUE)}
  	```


## 2. **using normalized and scaled assay for clustering and ploting**
DefaultyAssay(obj) <- "integrated" # for clustering and plotting

## 3. when need to subset or grouping, often need to **switch active object**
Indent(obj) <- "cell.type" # need if to define idents inside the FindMarker function


# about different types of assay
########## types of assays
1. the "RNA" assay contains the raw expression data in count slot and normalized log based in data slot, 
2. the "refAssay" is an optional reference dataset for integration,
3. the "integrated" assay holds the integrated expression data after the integration process; which only have ~ 2k HVG in its count/data matrix; only good for integration; not suitable for find markers and explore gene expression.
4. the "SCT" assay: a little bit tricky, a better alternative to std preprocessing wf: library size normalization, log transformaton, and scale; its counts in count slot has been done some correction

########## types of applications
## Any functions used for calculate or visualiza gene expression need to use RNA assay; default use the data slot or scale slot (DoHeatmap)
1. visualize gene expression; (violin plots, shading expression values on tSNE/UMAP. Anywhere you want to see the expression differences such as FeaturePlot, VinPlot), finding markers; calculate AVG....: need to use "RNA" assay; The data slot (normalized and log-based data) is the defaulty used.
2. DoHeatmap: an exception for item 1 abover; default to use scaled data (mean-centered, sd-adjusted) due to show the up and down using color code
## Any function related with dimension reduciton, not matter clustering or visualing the embedding, use "integrated" assay
3. umap, pac, Dimplot, etc visualization: use reducted dimension data coordinates, so can use "integrated" assay 
so if did integaration and/or SCT, need to go back RNA assay, perform LogNormalize and scale again then can precede to visulaizing gene expression and finding markers. default data will be normalized and log-based data, except Heatmap which need scaled data

Conclusion: 
1. use "RNA" for gene expression analysis, e.g., find marker genes (DE), visualize gene expression, AVG gene expression refer my sildes of "tips";
2. using normalized and scaled assay, e.g, "integrated" assay for clustering and reduced dimension plotting

so if use SCT doing normalization, still do not forget to perform LogNormalize (normalization and log transform) and scaling on the RNA assay (not on integrated assay), which will be used for visualizing gene expression and finding markers

# about integration
1. merge function is only applicable for technical replcates which do not have much batch effect
2. integration function is need for integrating datasets which have batch effect due to different reasons
3. Do not peform a second round of SCT on integrated object if objectes before integration have been subjected a first round of SCT
4. After integration, need to go back to RNA assay to perfrom a new round of normalization and scale again for the purpose of downstream calculation or visualizing gene expression.

# about different data slots
"RNA" assay
counts: Stores unnormalized data such as raw counts or TPMs
data: Normalized data matrix (default for FindMarker) # log based
scale.data: Scaled data matrix, these are basically just z-scores; not applicable for FindMarkers; BUT for DoHeatmap

"SCT" assay
counts: corrected counts
data: log1p(correced_counts)
scale.data: pearson residuals

So SCTransform's GetAssayData(obj, slot="data") is not equal to RNA's NormalizeData(obj).


# about seurat versions
To maintain compatibility with previous workflows, new Seurat objects will use the previous object structure by default
To use new Seurat v5 assays in older version seurat, e.g., v4.3 : Please run: options(Seurat.object.assay.version = 'v5') 
my current r4seurat conda env in biowulf already installed v5, so no need to define the options
Most functions are compatiable to v5, including Azimuth, SCTranform

# about seurat object: a s4 object
rows: gene
columns: cell
obj[genes, cells] #for slice a subset

# about seurat meta.data: a dataframe
rows: cell
columns: attributes