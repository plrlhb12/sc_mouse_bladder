# general info about seurat object

#####################################################
# S4 object and Seurat
# # Seurat is a S4 object containing slots which contains either a string, a factor, a s3 object, a dataframe (meta.data), a large list of S4 object, and etc.
#####################################################

#####################################################
# Assays, key objects in the Seurat object; one Seurat can store more than 1 assay, only one will be as default assay
# one assay is one slot for storing expression matrix
#####################################################


#####################################################
# Seurat v5 assays store **expression data** in layers in an assay slot
# layers: one assay is for storing expression matrix, it can have one or multiple layers
# it could have one layer of "counts" as the raw counts, OR be added with more layers when applied normalization and scale
# These layers can store as 
# 1. raw, un-normalized counts (layer='counts'), 
# 2. normalized data (layer='data'),
# 3. z-scored/variance-stabilized data (layer='scale.data')
#####################################################
Layers(hamsterlung[["RNA"]]) # extract the layer names from one assay

###################################################
# structure of the seurat object, a S4 object
# assays, graphs, reductions in a Seurat object are one or multiple S4 objects
# To access slots of an S4 object you use @, not $:
###################################################

# using pbmc3k as an example, which has two assays, but haven't subjected normalization, scale, and pca
pbmc3k@assays # A list of assays for this project
pbmc3k@meta.data # meta-information about each cell, starting with number of features detected (nFeature) and the original identity class (orig.ident); more information is added using AddMetaData
# Name of the active, or default, assay; settable using DefaultAssay
pbmc3k@active.assay
# The active cluster identity for this Seurat object; settable using Idents; default using the project name if not set 
pbmc3k@active.ident
# A list of Graph objects which contain S4 objects
pbmc3k@graphs
pbmc3k@neighbors
# A list of dimensional reduction S4 objects for this object, e.g., pca, umap
pbmc3k@reductions
pbmc3k@images
pbmc3k@project.name
pbmc3k@misc
pbmc3k@version
pbmc3k@commands
pbmc3k@tools
  
# using lung72h as an example
hamsterlung@assays
class(hamsterlung@assays)
colnames(hamsterlung@meta.data) # orig.ident" "nCount_RNA"  "nFeature_RNA"  "percent.mt" "RNA_snn_res.0.5" "seurat_clusters"
hamsterlung@active.assay
hamsterlung@active.ident #Levels: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
class(hamsterlung@active.ident) # "factor"


###################################################
# structure of the "graphs" slot
# a list of lists which contain s4 objects
# each sublist is a S4 object which can use @ to access slots, e.g, attributes
# # access each by either index (object@graphs[1]) or name (object[[Names]])
# what are i and p stand for??
###################################################

# dissect the first list (RNA_nn) inside the "graphs" slot
# access it by either index (object@graphs[1]) or name (object[[Names]])
length(hamsterlung@graphs) # a list of 2 lists; RNA_nn and RNA_snn
class(hamsterlung@graphs[1]) # a large "list" with 2 dimenions 11461x11461
graph1_lung <- hamsterlung@graphs[1]
class(graph1_lung) # a large list containing S4 objects
graph1_lung[["RNA_nn"]]@i[1:5] # [1]    0  107  258 1218 2785
graph1_lung[["RNA_nn"]]@p[1:5] # [1]  0 21 37 46 58

###################################################
# structure of the the "reduction" slot
# a list of lists
# each sublist is a S4 object which can use @ to access slots, e.g, attributes
# # access it by either index (object@graphs[1]) or name (object[[Names]])
# the default sequence of reduction:umap, tsne, then pca
###################################################

# dissect the first list (umap) inside the "reduction" slot
hamsterlung@reductions # a list of 2 lists (pca and umap)
pca_data <- hamsterlung@reductions[2]
pca_data
class(pca_data[["umap"]]@cell.embeddings) # "matrix" "array"
pca_data[["umap"]]@cell.embeddings
#                                   UMAP_1        UMAP_2
#  AAACAAGCAAGGTTATATGTTGAC-1 -4.715840162  2.9832330796


###################################################
# other slots of the S4 object hamsterlung
###################################################

hamsterlung@images
hamsterlung@project.name
hamsterlung@misc
hamsterlung@version
hamsterlung@commands # contain used commands
hamsterlung@tools


###################################################
# Method 1: access individual component stored in slots using the way of S4 object 
# using @ to access a slot AND then []/[[]] to slice a list or data.frame
###################################################

pbmc3k@assays[1] # get the 1st assay from the asssays slot
pbmc3k[["RNA"]] # alternative way

###################################################
# access the raw count matrix inside an assay: 
# rownames are genes; 
# each column is a  cell; it is Matrix of class "dgCMatrix"
###################################################
pbmc3k[["RNA"]]$counts[1:3, 1:3]
## 3 x 3 sparse Matrix of class "dgCMatrix"
##       AAACATACAACCAC AAACATTGAGCTAC AAACATTGATCAGC
## CD3D  4              .                         10 


###################################################
# access the metadata inside an assay: 
# meta.data: info associated with each cell
# rownames are cells; 
# each column is an atrribute about a cell;
###################################################
pbmc3k@meta.data[["orig.ident"]] # get the orig.indent from meta.data
pbmc3k@meta.data[1, 1:4]

##                   orig.ident nCount_RNA nFeature_RNA percent.mt
## AAACATACAACCAC-1     pbmc3k       2419          779  3.0177759

# Show QC metrics for the first 5 cells
head(pbmc3k@meta.data, 5)


###################################################
# Method 2: access individual component stored in slots using the way of Seurat Object which is based on S4 object 
# but has its own special functions
# get or store expression data, i.e., counts/matrix using object[[AssayName]]$counts/data/scaled.data
# get all metat.data useing object[[]]  
# get a single column from the meta.data using use object$colname 
###################################################

# get raw counts from the assay of "RNA"
pbmc3k[["RNA"]]$counts

# get normalized counts if there is one
pbmc3k[["RNA"]]$data

# change default assay
newAsssay <- subset(pbmc3k[["RNA"]], cells=columns(pbmc3k)[1:100])
DefaultAssay(pbmc3k) <- "newAssay"

# get meata data
pbmc3k[[]] # a dataframe
pbmc3k$seurat_annotations # extract a single features stored


############### get cell names/barcodes using Cells OR colnames
# only get those for the default assay
Cells(pbmc3k)
# get all for the entire object
colnames(pbmc3k)

############### sanity check the assay version or mem use or storage
class(pbmc3k[["RNA"]])
class(pbmc3k[["RNA"]]$counts)



