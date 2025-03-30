# here is my concise script for run annotation using Azimuth reference database in R
# rscript file.h5ad

################# if in bash
# input: either direct from the output dir from Cellranger or existing h5ad object
INPUT=${1}

if [[ "${INPUT##*.}" == "h5ad" ]]; then
	rscript $INPUT
elif [[ "${INPUT##*.}" == "h5" ]]; then
	#statements
elif [[ condition ]]; then
	#statements
else; then
	#statements
fi



################# if in r

######## read files and do mapping

library(Seurat)
library(SeuratData)
library(patchwork)
library(sctransform)
library(Azimuth)
#library(BPCells)
#library(dplyr)
#library(scater) # for more options of visualization
#library(ggplot2)
#library(SeuratDisk) # for interoperationality
options(Seurat.object.assay.version = "v5") # define using v5


######### load data and do reference mapping
file_name <- "your_file.h5ad"

if (grepl("\\.h5ad$", file_name)) {
  print("The file name ends with '.h5ad'")
  obj <- RunAzimuth(file_name, reference="lungref")
} else if (grepl("\\.h5$", file_name)) {
  print("The file name ends with '.h5'")
  obj <- Read10X_h5(file_name)
  obj <- RunAzimuth(obj, reference="lungref")
} else if (file.info(file_name)$isdir) {
  print("The file name is a directory")
  obj <- Read10X(file_name)
  obj <- RunAzimuth(obj, reference="lungref")
} else {
  print("The file name does not match any condition")
}


########## pipeline with SCTtransform 


