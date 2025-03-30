# convert count mtx matrix to 10x HDF5 file using the package of scCustomize
# the cons of this method: the gnerated h5 file will, the metadata's column would be a bit different form that from cellrange
library(scCustomize)
for (folder in count.folders){
  inpath <- paste0("../raw_data/count_folder/", folder)
  print(inpath)
  save.name <- "sample_filtered_feature_bc_matrix.h5"
  print(save.name)
  outpath <- paste0("../raw_data/h5files/", folder)
  print(outpath)
  dir.create(outpath)
  Create_10X_H5(raw_data_file_path = inpath, save_file_path = outpath, save_name = save.name)
  #becasue the generated h5 has extra characters as sample_filtered_feature_bc_matrix.h5_1e0588ad457d9.h5, rename them
  h5name <- list.files(outpath, pattern = ".h5$")
  print(paste0("h5file is named as", h5name))
  oldname <- paste0(outpath, "/", h5name)
  print(oldname)
  newname <- gsub("\\.h5_.*h5$", ".h5", oldname)
  print(newname)
  print(paste0("rename as ", newname))
  file.rename(from = oldname, to = newname)
}