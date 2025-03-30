# The integrated obj takes 20G mem, and sinteractive --mem=100g is still challenge; so increase to 240g; increse the max mem in .Renviron up to 200g

# sinteractive --mem=240g --cpus-per-task=12 --gres=lsratch:100
# cd $IRF
# source myconda
# conda active r4seurat
# module load rstudio R
# rstudio&
# setwd("/vf/users/pengl7/IRF/seurat/flex2023/")
# mkdir "20240119-integrate"
# create a new project using the folder of 20240119-integrate

library(Seurat)
library(SeuratData)
library(patchwork)
library(Azimuth)
library(sctransform)
library(ggplot2)
library(hdf5r)
library(dplyr)
library(tidyr)
library(DESeq2) # if want to perform DE using this option
library(reticulate)
use_condaenv("/data/pengl7/conda/envs/r4seurat")
py_config()
# need to first create a conda env as r4seurat; conda active r4seurat; conda install leidenalg; conda install numpy pandas; start rstudio&
leidenalg <- import("leidenalg")