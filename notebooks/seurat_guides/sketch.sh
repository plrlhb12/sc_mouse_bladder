# https://satijalab.org/seurat/articles/seurat5_sketch_analysis
# Sketch-based analysis in Seurat v5
Infrastructure for on-disk storage of large single-cell datasets instead of storing expression matrices in memory
‘Sketching’ methods to subsample cells from large datasets while preserving rare populations
after sketching, the subsampled cells can be stored in-memory, allowing for interactive and rapid visualization and exploration.
We store sketched cells (in-memory) and the full dataset (on-disk) as two assays in the same Seurat object. Users can then easily switch between the two versions, 
providing the flexibility to perform quick analyses on a subset of cells in-memory, while retaining access to the full dataset on-disk.
ultilize BPCells  utilizes bit-packing compression to store count matices on-disk, and optimized streaming-compatible C++ code to substantially improve I/O and computational performance when working with on-disk data. 

library(Seurat)
library(BPCells)
library(ggplot2)
# needs to be set for large dataset analysis
options(future.globals.maxSize = 1e9)

1. create this on-disk representation using BPCells: Remember that one matric of h5 already has row and column names, no more other info. 
# Read a feature matrices from 10x
if start from a h5 file, open_matrix_10x_hdf5 can be used to read it 
	brain.data <- open_matrix_10x_hdf5(path = "/brahms/hartmana/vignette_data/1M_neurons_filtered_gene_bc_matrices_h5.h5")
	# Write the matrix to a directory
 	write_matrix_dir(
   		mat = brain.data,
   		dir = '/brahms/hartmana/vignette_data/bpcells/brain_counts')
	# load the matrix from disk (this only generates the connection, doesn't load the whole data) using open_matrix_dir
	brain.mat <- open_matrix_dir(dir = "/brahms/hartmana/vignette_data/bpcells/brain_counts")
	brain.mat <- Azimuth:::ConvertEnsembleToSymbol(mat = brain.mat, species = "mouse")
	brain <- CreateSeuratObject(counts = brain.mat) # will this loose all the meta info??

if already have a seurat object: convert the counts matrix of your RNA assay to a BPCells matrix
	# # Write the counts layer to a directory
 	write_matrix_dir(mat = obj[["RNA"]]$counts, dir = '/brahms/hartmana/vignette_data/bpcells/brain_counts')
 	counts.mat <- open_matrix_dir(dir = "/brahms/hartmana/vignette_data/bpcells/brain_counts")
	obj[["RNA"]]$counts <- counts.mat #replace the orignial count layer in-mem with the BPCells count matrix on-disk