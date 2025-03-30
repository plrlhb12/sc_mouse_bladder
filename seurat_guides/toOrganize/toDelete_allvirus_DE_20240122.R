
# peform DE among days by cell types
# using For loop for different cell types
```{r}
cell.types <- levels(factor(allvirus$cell.type)) 
days <- c("D03", "D07", "D10")
Idents(allvirus) <- "cell.type" 

for (type in cell.types){
  sub.obj <- subset(allvirus, idents = type)
  run.name <- type
  dir.create(paste0("outputs/",run.name))
  Idents(sub.obj) <- sub.obj$days
  # Assuming you have a list of 7 objects named 'compare.list'
  compare.list <- days
  num_objects <- length(compare.list)
  
  # Loop through each pair of objects and apply your comparison function
  for (i in 1:(num_objects - 1)) {
    for (j in (i + 1):num_objects) {
      # Skip self-comparisons
      if (i != j) {
        compare <- paste0(compare.list[[i]], "_vs_", compare.list[[j]])
        ident1 = compare.list[[i]]
        ident2 = compare.list[[j]]
        markers <- FindMarkers(sub.obj, ident.1 = ident1, ident.2 = ident2)
        write.csv(markers, file = paste0("outputs/", run.name, "/", run.name,"_", compare,"_DE_genes.csv"))
        
        top50 <- markers %>%
          slice_max(n = 50, order_by = avg_log2FC)
        write.csv(top50, file = paste0("outputs/", run.name, "/", run.name, "_", compare,"_top50_genes.csv"))
        
        bottom50 <- markers %>%
          slice_min(n = 50, order_by = avg_log2FC)
        write.csv(bottom50, file = paste0("outputs/", run.name, "/", run.name,"_", compare,"_bottom50_genes.csv"))
      }
    }
  }
  
  # continue to compare D3 to the others
  compare <- "D3vsD7andD10"
  markers <- FindMarkers(sub.obj, ident.1 = "D03")
  write.csv(markers, file = paste0("outputs/", run.name, "/", run.name,"_", compare,"_DE_genes.csv"))
  
  top50 <- markers %>%
    slice_max(n = 50, order_by = avg_log2FC)
  write.csv(top50, file = paste0("outputs/", run.name, "/", run.name, "_", compare,"_top50_genes.csv"))
  
  bottom50 <- markers %>%
    slice_min(n = 50, order_by = avg_log2FC)
  write.csv(bottom50, file = paste0("outputs/", run.name, "/", run.name,"_", compare,"_bottom50_genes.csv"))
}
```
# wrap the above block as a function so I can call it frequently
# perform_DE_compareDay_amongCellType
# have done on allvirus, need to do on allctrl # 20240122
```{r}
# paramters: a seurat object, a folder to becreated; example c(allvirus, allvirus)
# make sure 
perform_DE_compareDays_amongCelltype <- function(allvirus, outfolder) {
  wholedata <- allvirus
  # Set cell types
  days <- levels(factor(wholedata$days))
  cell.types <- levels(factor(wholedata$cell.type))
  # make sure create the new output folder "outputs/newOutFolder"
  parent.dir <- paste0("outputs/", outfolder)
  if (!dir.exists("outputs/")){
    dir.create("outputs/")
    print("generate the folder of outputs/")
  }
  dir.create(parent.dir)
  print(paste0("generate the folder of ", parent.dir))
  
  # looping all cell types
  Idents(wholedata) <- "cell.type"
  for (type in cell.types){
    sub.obj <- subset(wholedata, idents = type)
    run.name <- type
    subfolder <- paste0(parent.dir, "/", run.name)
    dir.create(subfolder)
    
    # compare days
    Idents(sub.obj) <- sub.obj$days
    compare_list <- days
    num_objects <- length(compare_list)
    for (i in 1:(num_objects - 1)) {
      for (j in (i + 1):num_objects) {
        if (i != j) {
          compare <- paste0(compare_list[[i]], "vs", compare_list[[j]])
          ident1 = compare_list[[i]]
          ident2 = compare_list[[j]]
          markers <- FindMarkers(sub.obj, ident.1 = ident1, ident.2 = ident2)
          write.csv(markers, file = paste0(subfolder, "/", run.name, "_", compare, "_DE_genes.csv"))
          
          top50 <- markers %>%
            slice_max(n = 50, order_by = avg_log2FC)
          write.csv(top50, file = paste0(subfolder, "/", run.name, "_", compare, "_top50_genes.csv"))
          
          bottom50 <- markers %>%
            slice_min(n = 50, order_by = avg_log2FC)
          write.csv(bottom50, file = paste0(subfolder, "/", run.name, "_", compare, "_bottom50_genes.csv"))
        }
      }
    }
    
    # Continue to compare D3 to the others
    compare <- "D3vsD7andD10"
    markers <- FindMarkers(sub.obj, ident.1 = "D03")
    write.csv(markers, file = paste0(subfolder, "/", run.name, "_", compare, "_DE_genes.csv"))
    
    top50 <- markers %>%
      slice_max(n = 50, order_by = avg_log2FC)
    write.csv(top50, file = paste0(subfolder, "/", run.name, "_", compare, "_top50_genes.csv"))
    
    bottom50 <- markers %>%
      slice_min(n = 50, order_by = avg_log2FC)
    write.csv(bottom50, file = paste0(subfolder, "/", run.name, "_", compare, "_bottom50_genes.csv"))
  }
}
```