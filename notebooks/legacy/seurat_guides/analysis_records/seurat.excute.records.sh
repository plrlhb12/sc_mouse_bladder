
####################################### performed sct and annotate for indiviual project (38 in total) in 2023Aug using all individual h5 file
refer to cellranger_record.sh

Note that: 
# 29Aug2023 found 10 h5 files from pool789 look suspicious (some files' sizes are very small; some of them were revealed to haven't been updated with the output generated on 
sars-gene-incorportaed barcode file)
# recopy the 10 h5 from original cellranger output folder in Sep6 2023, but accidently deleted all h5 in the origianl folder
# copy the original barcode files to Biowulf so that I can regenerate h5 file in the future (didn't copy the 2 from pool9, didn't recall why; in my impression, the h5 of them didn't give any trouble before)

####################################### 20230911-12
For the allvirus object:
integrate from all individual rds which has been SCT and annotated before

starting from RNA assay
folowing Seurat's integration for large dataset cPCA instead of cca
re-annotate using integrated assay
however, didn't peform cluster (error out)

# issues encounterd during seuarat 20230831
# cann't perform clustering for the integrated object 
# allocating vector of size 230.2
# Error: vector memory exhausted (limit reached?)
#Error in unlist(object) :
#long vectors not supported yet: ../../src/include/Rinlinedfuns.h:537

#####################################################################
# 20230119: organize and continue sc anaysis
1. comb what I have done from Aug to Sep 2023, re-organize folders, add timeline.exls
2. modify integrate_allvirus.Rmd of 20230912 to remove those scripts which were executed before 20230911 and only kept the commands exectuted only in 20230912.
3. copy and rename the script as integrate_virus_20230119.Rmd using the existing rds (allvirus_az.rds) to generate count by cell type by day. have not remove unrelavent commands.

# 20240120: since I have to integrate all_ctrl and all_virus together, I decide to start from begnining so I can re-warm what I have done and remove the meta and parts not requred and relavant
1. using the earlies script of integrate.Rmd of 20230906, rename as 20230120_all.Rmd
sinteractive --mem=240g --cpus-per-task=16 --gres=lscratch:100
source conda
conda activate r4seuratset
module load rstudio R
rstudio&
getwd()
setwd("/vf/users/pengl7/IRF")
#click go to workdirectory in rstudio
cd seruat/flex2023
#create a new R project 20230120_all in rstudio, chose create a git repo
cp ../20230831-integrate/integrate.Rmd ./20230120_all.Rmd # execute in terminal


# 20240122: modify DE assay; save them into seperate scripts for future call
# risk: not sure change the defaulty assay of the allvirus or not

source("perform_DE.R") # the funcitons are valid, but the way to call it haven't been tested
peform_DE(allvirus, "allvirus")

#20240123: peform DE assay for allctrl.RDS; modify and add metadata; saved the new rds and new script (integrate_ctrl_20240123.Rmd)
#20240123-24: peformed DE analysis using integrate_virus_20240119_good.Rmd for allvirus.RDS using the above source("peform_DE") ways to redo and newdo all the DE tasks as listed in task.ppt
saved some functions in seperate R scripts for future faster reference; rerun allvirus samples's DE due to DefaultyAssay issue; saved old output as DE_notRNAsasay
(# may do for allvirus.rds: add one more metadata column: predicted.cell.type which was prevous changed to cell.type)
#20240124:tested integrate allctrl and allvirus but fail as in the folder of integrate_all_20240120_failed; submit a issue report to seurat git repo
2. encounter error in merge and integrate
https://github.com/satijalab/seurat/issues/8369

####################### 20240125
# previous said 10 problemtic h5 were actually deleted before Aug6 2023. All the h5 inisde raw_data/h5files generated on Sep6 2023, include the 8 in the folder of "wrong" should be all correct
# 20240125: try to re-generate h5 file from count folder using scCustomize in r4seurat conda enviroment; notice that here .hs means 5HDF
# found only 8 projects's count folders intead of 10 have been transfered into biwoulf. 
# The other 2 projects (4M-Omicron-D7 and D10) .h5 file are inside the folder of h5file and not in "wrong" folder, but their original h5 in Cellrange output have been deleted. 
# For the sake of safety, I updated all 10 from pool7,8,9 by regenerate them from counts (they are a little bit larger than the old ones). and keep their older copies in a seperate folder
# copy paste all 10 h5 files back to their original Cellrange output folder in Locus.

################## 20240126: rename h5 files; but accidently cause some dataset overwirite each other ; recopy from snapshot yesterday to here
re-read all h5 files and only do minimal preprocessing and save them. 2 WA1-D03 AND D07 samples casue issue in opening the H5; 
recopy h5 from locus pool 10; found the folder name haven't corrected from 4D-Omicron-D3/7 to 4D-WA1-D3/7; correced the folder names in pool 10
Try to merge all object as one


########################### 20240129: merged all objectes as 1 using scCustomize; fail during IntegrateLayers OR findClusters, highest mem use was 412 G, exceed the Renvion defined 400G


#################### 20240124 update on integration: try to integrate allCtrl and allVirus together instead starting from begining
4. After integration, need to go back to RNA assay to perfrom a new round of normalization and scale again for the purpose of downstream calculation or visualizing gene expression.


#################### 20240129 finially integrated all togther, which haven't gone through clustering and annotation
individual h5 - normalization - rds - one copy of merged.rds in integrated20240126/rds
individual h5 - sct --- rds - one copy of integrated.rds -- save in integrate.rds in integrated20240130
for (i in 1:length(ids)){
  cat("processing ", ids[i])
  obj <- readRDS(rds.path[i]) %>%
    SCTransform() %>%
    RunPCA()%>%
    saveRDS(file = paste0(parent.path, "/", ids[i], ".rds"))
  }

# perform integraton 
features <- SelectIntegrationFeatures(object.list = obj.list)
anchors <- FindIntegrationAnchors(object.list = obj.list, reference = c(19, 24, 35), reduction = "rpca", dims = 1:50)
integrated <- IntegrateData(anchorset = anchors, dims = 1:50, normalization.method = "SCT")
saveRDS(integrated, file = paste0(parent.path, "/", "integrated.rds"))

#################### 20240130 generated the on-disk of bpcells, but the file size didn't reduced in both rstuido mem use and on disk usage
save a copy of on-disk integrated.rds inside 20240130_skecth

#################### 20240131 try to use sketch strategy 
1. On tutorial dataset: no problem
2. on my own integrated.rds on-disk version, failed on 2 steps: 
        1. RunAzimuth; skip this step
        2. ProjectData; still failed as below


obj <- readRDS("integrated.rds")
format(object.size(obj), units = "Gb")

DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- SketchData(
  object = obj,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

print(obj)

DefaultAssay(obj) <- "sketch"
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:50)
obj <- FindClusters(obj, selection.method = "vst", nfeatures = 2000, algorithm = 4)
obj <- RunUMAP(obj, dims = 1:50, return.model = T)


obj <- ProjectData(
  object = obj,
  assay = "RNA",
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  umap.model = "umap",
  dims = 1:50,
  refdata = list(cluster_full = "seurat_clusters")
)

Error in dimnames(x) <- dn : 
  length of 'dimnames' [2] not equal to array extent
In addition: Warning message:
In matrix(data = unlist(proj.pca.list), nrow = nrow(proj.pca.list[[1]]),  :
  data length [455900] is not a sub-multiple or multiple of the number of columns [222935]


# now that integration is complete, rejoin layers
obj[["RNA"]] <- JoinLayers(obj)

######## back to orignal method of clustering and annotation
DefaultAssay(integrated) <- "RNA"
integrated[["RNA"]] <- JoinLayers(integrated, assay = "RNA") 
Error in `[<-.data.frame`(`*tmp*`, , i, value = new("Seurat", assays = list( : 
  replacement has 12943 rows, data has 222935


################### 
problem in reading 4D-Ctrl-D7.h5 ??? 


###### use top to check the cpus and mem use in terminal

#################### 20240202_sketch
20240202_sketch.Rmd
20240203_sketch_cluster.Rmd------ good Rmd comments
20240206_sketch_cluster.Rmd------ good Rmd comments

merged.rds  object just created from bpcell matrix, haven't been sketched yet
sketch.rds  object just finished sketch, normal, findingvaribles post-sketch,  but not integrating layers yet
integrated.rds  object finished integrating layers and normalized but not finding clusters
integrated_llouvein.rds object finished louvein clustering
integrated_projected.rds  

summary the error I encounted in prevously numerious trying
1. orignal 38 h5 files will give different formats on bpcell data list: the original 28 from cellranger use ensemble id as feature ID; 
  while the 10 h5 converted from scCustmize pakage use gene symbol as feature ID;
  This discrepency cause two sets features, each has 19065 genes; thus there isn't any variable features found duing the sketch workflow on integration
  SOLUTIONS: regenerate h5 files using cellranger

2. cell barcode names are not unique: error occured during sketch integrateLayer(): becasue more than 2 sampelesets use the same cell barcodes, which causet uncompatiable
  number of dimentsion.
  sulutions: add suffix to cellbarcode names in each bpcell matrix

2. long vector using leiden clustering: 
  solutions: change method to "igraph"

3. azimuth problem: ! GetAssayData doesn't work for multiple layers in v5 assay.
  solutions: only perform it after joinLayers()


####################### 20240206_newinstall
#since sctransform and azimuth didn't work as expected and doubt the verison of seuart; tried to install all related packages;
only succeeded in 3 as below; but failed in install azimuth, bpcells and signac due to git credential; the uncompletted intallation
may affecte the apropriate use of the pre-installed packages; have to deleted those starig with ""00LOCK-"".

remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
remotes::install_github("satijalab/seurat-data", "seurat5", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", "seurat5", quiet = TRUE)

#################### 20240131_integrate_oldway (executed on 20240206)
#lots of try and error on Azimuth and sctransform;
#remove and reinstall some packages;
#SCTransform failed to work after runing once;
#can ignore this folder;20240206_integrate.Rmd
#lots of try and error on Azimuth and sctransform;remove and reinstall some packages;
#SCTransform failed to work after runing once;can ignore this folder;
#20240206_integrate.Rmd

####################### 20240210_repeatindividuals-------partially successful
#20240210: only sucessfuly runing std normalize and annotate all individuals h5 using std normalization methods (failed in using sctrasnfrorm)
#Haven't done any integration using this batch of rds files
#rds (individuals-20240210) only std normaliza and annotated using new h5 generated on 20240202


###################### 20240212_integrateAll, maxmem ~ 200G --------- succesffuly
#succesffuly interate allctrl and allvirus using the old integration method (not working using new integration IntegrateLayer())
#"change the name of integrated_dr" which was generated during RunAzimuth() process using reference mapping or delete it make this successful
#save as integrated_all.rds in which meta.data has lots of columns that control doens't have values`


###################### 20240220 contiune to the project of above "20240212_integrateAll" 
# perform the left analysis; maxmem 50G ------ sucesful
# deleted several columns of metadata and generated a new column of "groups" saved as **all_cleaned.rds** from **integrated_all.rds**
# calculate DE and AVG
# generated skyline endpoint and transfered this output to my folder in /data/irf/gn/Lirong/flex_20240220/all


################### 20240328-repeateColton try to copy Colton's solutions
# failed to install package in Biowulf in the new conda env
# created new env and package in skyline but failed to appy Aziumth and fail to connect it to Rstudio


################### 20240329-subset: 
# add a little more stringent filtering on the finial rds (# SUBSET by only keep cells which have more than 3000 gene and less than 300; 
# the prevous filtering threthold is low (min.cells = 3, min.features = 200)) from the all_cleaned.rds to 20240329_filtered.rds in 20240329 folder
# tried trajectory
# tried doublets calculation
# sub-cluster unclear cell types

################### 20240424-re-DEG
# regenerate DEG list for pathway analysis
# follow the method used in the scripts of 20240220_integarteAll_continueDEavg in 20240212-integrateAll_ContrieDE_good using separate functional scripts
# since all my previous cell number and graph were generated on all_cleaned.rds (i.e., integarted_all.rds). so keep using it instead of the one 20240329_filtered.rds





