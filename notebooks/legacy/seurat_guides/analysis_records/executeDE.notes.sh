# perform DE by days
# It usually recommend using "RNA"assay, although, Seurat also say SCT could do this purpose

#################################################
# most important tips in DE assay and AVG gene expression assay
##################################################

1. Make sure set the active assy is RNA; best to define it outside a function if want the change be permenant
* please note that even although I change the defaulty assay to RNA inside the function of perform_DE, 
* the change only occured inside the function; it will not return the modified object; so outside the funciton, the obj still have the previous default assay

2. make sure reset active.idents before execution

################################################## 20240119-24
# I have generated separte scripts for the functions which using for loops to DE on conditions
1. pefrome one comparsion or several comparisons on 1 condition; saved as perform_DE_***.Rmd


################################################## history
# perform DE by days on integrate_virus.Rmd on Sep2023
It usually recommend using "RNA"assay, However, Seurat also say SCT could do this purpose
# perform DE using SCT data
# MM is interested in export 100 DE for each cell types
# here I changed the seurat_clusters to the value of celltype, save a copy of it at the same time

# Perform DE analysis using SCT # not sure whether I used it successfuly or not
https://satijalab.org/seurat/reference/prepsctfindmarkers
NEED PrepSCTFindMarkers()
Idents(sub.obj) <- sub.obj$days
sub.obj <- PrepSCTFindMarkers(object = sub.obj) # takes about 1 h
