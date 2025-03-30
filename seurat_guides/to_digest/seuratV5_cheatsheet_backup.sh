# https://satijalab.org/seurat/articles/seurat5_essential_commands.html
library(Seurat)
library(SeuratData)
library(BPCells)
library(dplyr)
options(Seurat.object.assay.version = "v5") # define using v5

pbmc3k <- LoadData("pbmc3k") # from build-in dataset
mousebrain1m <- readRDS("/brahms/hartmana/vignette_data/1p3_million_mouse_brain.rds") # from seurat object RSD

#  RNA assay is of the Assay5 class
class(pbmc3k[["RNA"]])
class(mousebrain1m[["RNA"]])

###########  Access and store **expression data** using The $ and double-bracket [[]] symbols
# access the counts matrix from the RNA assay
counts_matrix <- pbmc3k[["RNA"]]$counts

# Add a layer Equivalent to running pbmc3k <-NormalizeData(pbmc3k)
pbmc3k[["RNA"]]$data <- NormalizeData(pbmc3k[["RNA"]]$counts)

# Delete a layer
pbmc3k[["RNA"]]$data <- NULL

# pbmc3k counts matrix is stored in-memory as shown "BPCells"
class(pbmc3k[["RNA"]]$counts)


########## Access cell names using Cells() or colnames() and metadata using [[]] and $
pbmc3k[["RNAsub"]] <- subset(pbmc3k[["RNA"]], cells = colnames(pbmc3k)[1:100])
DefaultAssay(pbmc3k) <- "RNAsub"
length(Cells(pbmc3k))

# get all object metadata
pbmc_metadata <- pbmc3k[[]]
# get list of metadata columns
colnames(pbmc_metadata)
# get annotations stored in metadata
annotations <- pbmc3k$seurat_annotations


############ Create Seurat or Assay objects; can define version
# create v3 object
options(Seurat.object.assay.version = "v3")
pbmc.counts <- Read10X(data.dir = "/brahms/hartmana/vignette_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.counts)
class(pbmc[["RNA"]])
# create v5 object
options(Seurat.object.assay.version = "v5")
pbmc.counts <- Read10X(data.dir = "/brahms/hartmana/vignette_data/pbmc3k/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.counts)
class(pbmc[["RNA"]])

# create a v3 assay
assay.v3 <- CreateAssayObject(counts = pbmc.counts)
# create a v5 assay
assay.v5 <- CreateAssay5Object(counts = pbmc.counts)
class(assay.v3)
class(assay.v5)

# create an assay using only normalized data and generate a object on it
assay.v5 <- CreateAssay5Object(data = log1p(pbmc.counts))
# create a Seurat object based on this assay
pbmc3k_slim <- CreateSeuratObject(assay.v5)
pbmc3k_slim

# convert a v5 assay to a v3 assay
pbmc3k[["RNA3"]] <- as(object = pbmc3k[["RNA"]], Class = "Assay")
# convert a v3 assay to a v5 assay
pbmc3k[["RNA5"]] <- as(object = pbmc3k[["RNA3"]], Class = "Assay5")


############# Working with layers, These layers can store 
# 1. raw, un-normalized counts (layer='counts'), 
# 2. normalized data (layer='data'), or 
# 3. z-scored/variance-stabilized data (layer='scale.data').
# by default, creates an RNA assay with a counts layer
obj <- CreateSeuratObject(counts = pbmc.counts)
# creates a normalized data layer
obj <- NormalizeData(obj, verbose = FALSE)
# extract only the layer names from an assay
Layers(obj[["RNA"]])

# split the layers into groups
................................


############### Accessing additional data
pbmc3k <- FindVariableFeatures(pbmc3k, verbose = FALSE)
pbmc3k <- ScaleData(pbmc3k, verbose = FALSE)
pbmc3k <- RunPCA(pbmc3k, verbose = FALSE)

# return variable features
# returns information from both assay (here from "count" layer), features (MS4A from the assay of "rna_"), cell embeddings ("PC_1") and meta.data ("nCount_RNA") as a data.frame
fetch_df <- FetchData(object = pbmc3k, layer = "counts", vars = c("rna_MS4A1", "PC_1", "nCount_RNA"))
head(fetch_df)

# get cell embeddings
head(Embeddings(object = pbmc3k[["pca"]])[, 1:5])

# get feature loadings
head(Loadings(object = pbmc3k[["pca"]])[, 1:5])



################## ################## ################## ##################
# New funcitons coming with V5
# 1. supervised mapping for annotation


################## ################## ################## ##################





################## ################## ################## ##################
# New funcitons coming with V5
# 2. sketch integration


################## ################## ################## ##################








################## ################## ################## ##################
# New funcitons coming with V5
# 3. Process large data using BPCells which work with Seurat Oject in mem while accessing counts on disk
################## ################## ################## ##################
devtools::install_github("bnprks/BPCells")
library(BPCells)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(Azimuth)
options(Seurat.object.assay.version = "v5")

# example 1: read from h5
brain.data <- open_matrix_10x_hdf5(path = "/brahms/hartmana/vignette_data/1M_neurons_filtered_gene_bc_matrices_h5.h5")
# Write the matrix to a directory, brain.data becomes an on-disk layer instead of in-memory layer saved in ..../brain_counts
write_matrix_dir(
  mat = brain.data,
  dir = "/brahms/hartmana/vignette_data/bpcells/brain_counts"
)
# Now that we have the matrix on disk, we can load it from the above saved on-disk layer instead of orignail h5 directory
brain.mat <- open_matrix_dir(dir = "/brahms/hartmana/vignette_data/bpcells/brain_counts")
brain.mat <- Azimuth:::ConvertEnsembleToSymbol(mat = brain.mat, species = "mouse")
brain <- CreateSeuratObject(counts = brain.mat)

# example 2: if alread have a Seurat object
obj <- readRDS("/path/to/reference.rds")
# Write the counts layer to a directory
write_matrix_dir(mat = obj[["RNA"]]$counts, dir = "/brahms/hartmana/vignette_data/bpcells/brain_counts")
counts.mat <- open_matrix_dir(dir = "/brahms/hartmana/vignette_data/bpcells/brain_counts")
obj[["RNA"]]$counts <- counts.mat

# if needed can convert a layer an on-disk matrix to in-memory matrix using as()
brain <- subset(brain, downsample = 1000)
brain[["RNA"]]$counts <- as(object = brain[["RNA"]]$counts, Class = "dgCMatrix")

# Saving Seurat objects in the dir where the on-disk layer was saved in "/brahms/hartmana/vignette_data/bpcells/" as previous for future loading, 
# brain.data in "/brahms/hartmana/vignette_data/bpcells/brain_counts/
# obj.Rds in "/brahms/hartmana/vignette_data/bpcells/brain_object/
saveRDS(
  object = brain,
  file = "obj.Rds",
  destdir = "/brahms/hartmana/vignette_data/bpcells/brain_object"
)


########################## example 3.  Load data from multiple h5ad files
file.dir <- "/brahms/hartmana/vignette_data/h5ad_files/"
files.set <- c("ahern_pbmc.h5ad", "jin_pbmc.h5ad", "yoshida_pbmc.h5ad")

# Loop through h5ad files and output BPCells matrices on-disk
data.list <- c()
metadata.list <- c()

for (i in 1:length(files.set)) {
  path <- paste0(file.dir, files.set[i])
  data <- open_matrix_anndata_hdf5(path)
  write_matrix_dir(
    mat = data,
    dir = paste0(gsub(".h5ad", "", path), "_BP")
  )
  # Load in BP matrices
  mat <- open_matrix_dir(dir = paste0(gsub(".h5ad", "", path), "_BP"))
  mat <- Azimuth:::ConvertEnsembleToSymbol(mat = mat, species = "human")
  # Get metadata
  metadata.list[[i]] <- LoadH5ADobs(path = path)
  data.list[[i]] <- mat
}
# Name layers
names(data.list) <- c("ahern", "jin", "yoshida")

# Add Metadata
for (i in 1:length(metadata.list)) {
  metadata.list[[i]]$publication <- names(data.list)[i]
}
metadata.list <- lapply(metadata.list, function(x) {
  x <- x[, c("publication", "sex", "cell_type", "donor_id", "disease")]
  return(x)
})
metadata <- Reduce(rbind, metadata.list)

options(Seurat.object.assay.version = "v5")
merged.object <- CreateSeuratObject(counts = data.list, meta.data = metadata)

saveRDS(
  object = merged.object,
  file = "obj.Rds",
  destdir = "/brahms/hartmana/vignette_data/bpcells/merged_object"
)


#################################### example 4: load big h5ad from Parse Biosciences
parse.data <- open_matrix_anndata_hdf5("/brahms/hartmana/vignette_data/h5ad_files/ParseBio_PBMC.h5ad")
write_matrix_dir(mat = parse.data, dir = "/brahms/hartmana/vignette_data/bpcells/parse_1m_pbmc")
parse.mat <- open_matrix_dir(dir = "/brahms/hartmana/vignette_data/bpcells/parse_1m_pbmc")
metadata <- readRDS("/brahms/hartmana/vignette_data/ParseBio_PBMC_meta.rds")
metadata$disease <- sapply(strsplit(x = metadata$sample, split = "_"), "[", 1)
parse.object <- CreateSeuratObject(counts = parse.mat, meta.data = metadata)
saveRDS(
  object = parse.object,
  file = "obj.Rds",
  destdir = "/brahms/hartmana/vignette_data/bpcells/parse_object"
)