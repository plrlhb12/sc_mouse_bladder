# paramters: a seurat object, a folder to create; example c(allvirus, allvirus)
# compare days for each celltype
perform_DE_compareDays_amongCelltypes <- function(allvirus, outfolder) {
  # makesure use RNA assay instead of integrated which may be the default
  wholedata <- allvirus 
  DefaultAssay(wholedata) <- "RNA"
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
          markers <- FindMarkers(sub.obj, ident.1 = ident1, ident.2 = ident2, logfc.threshold = 0, min.pct = 0)
          write.csv(markers, file = paste0(subfolder, "/", run.name, "_", compare, "_DE_genes.csv"))
          
          # top50 <- markers %>%
          #   slice_max(n = 50, order_by = avg_log2FC)
          # write.csv(top50, file = paste0(subfolder, "/", run.name, "_", compare, "_top50_genes.csv"))
          
          # bottom50 <- markers %>%
          #   slice_min(n = 50, order_by = avg_log2FC)
          # write.csv(bottom50, file = paste0(subfolder, "/", run.name, "_", compare, "_bottom50_genes.csv"))
        }
      }
    }
    
    # Continue to compare D3 to the others
    compare <- "D3vsD7andD10"
    markers <- FindMarkers(sub.obj, ident.1 = "D03", logfc.threshold = 0, min.pct = 0)
    write.csv(markers, file = paste0(subfolder, "/", run.name, "_", compare, "_DE_genes.csv"))
    
    # top50 <- markers %>%
    #   slice_max(n = 50, order_by = avg_log2FC)
    # write.csv(top50, file = paste0(subfolder, "/", run.name, "_", compare, "_top50_genes.csv"))
    
    # bottom50 <- markers %>%
    #   slice_min(n = 50, order_by = avg_log2FC)
    # write.csv(bottom50, file = paste0(subfolder, "/", run.name, "_", compare, "_bottom50_genes.csv"))
  }
}

# perform_DE_compareDays_amongCelltype(allctrl, "allctrl")
# perform_DE_compareDays_amongCelltype(allvirus, "allvirus")