# perform DE but compare between days for each cell type of each strain; 
perform_DE_compareDays_onCelltypeStrain <- function(allvirus, outfolder) {
  # set variables
  whole.obj <- allvirus
  upfolder <- outfolder
  print(upfolder)
  DefaultAssay(whole.obj) <- "RNA" # very important
  cell.types <- levels(factor(whole.obj$cell.type))
  days <- levels(factor(whole.obj$days))
  strains <- levels(factor(whole.obj$strains))

  # start strain loop
  Idents(whole.obj) <- "strains"
  for (strain in strains){
    strain.obj <- subset(whole.obj, idents = strain)
    parent_dir <- paste0(upfolder, "/", "byStrainCelltype")# here may need to customize
    if (!dir.exists(parent_dir)){
      dir.create(parent_dir)
    }
    dir.create(paste0(parent_dir, "/", strain)) # here may need to customize
    tryCatch({
      # start cell type loop
      Idents(strain.obj) <- "cell.type"
      for (type in cell.types){
        print(paste0("looping ", type, "for ", strain))
        # add error handle
        tryCatch({
          sub.obj <- subset(strain.obj, idents = type)
          run.name <- type
          final.out.dir <- paste0(parent_dir, "/", strain, "/", run.name)
          dir.create(final.out.dir)
          
          # set compare target of interest, here is to compare days
          compare_list <- days
          print(paste0("compare: ", compare_list))
          num_objects <- length(compare_list)
          
          Idents(sub.obj) <- "days"
          for (i in 1:(num_objects - 1)) {
            for (j in (i + 1):num_objects) {
              if (i != j) {
                compare <- paste0(compare_list[[i]], "vs", compare_list[[j]])
                print(compare)
                ident1 = compare_list[[i]]
                ident2 = compare_list[[j]]
                
                # Use tryCatch to handle the error
                tryCatch({
                  markers <- FindMarkers(sub.obj, ident.1 = ident1, ident.2 = ident2, logfc.threshold = 0, min.pct = 0)
                  write.csv(markers, file = paste0(final.out.dir, "/", run.name, "_", compare, "_DE_genes.csv"))
                  
                  # top50 <- markers %>%
                  #   slice_max(n = 50, order_by = avg_log2FC)
                  # write.csv(top50, file = paste0(final.out.dir, "/", run.name, "_", compare, "_top50_genes.csv"))
                  
                  # bottom50 <- markers %>%
                  #   slice_min(n = 50, order_by = avg_log2FC)
                  # write.csv(bottom50, file = paste0(final.out.dir, "/", run.name, "_", compare, "_bottom50_genes.csv"))
                  
                }, error = function(e) {
                  # Print a message or take other actions upon error
                  cat("Error:", e$message, "\n")
                })
              }
            }
          }
          
        }, error = function(e) {
          # Print a message or take other actions upon error
          cat("Error:", e$message, "\n")
        })
      }
      
    }, error = function(e) {
      # Print a message or take other actions upon error
      cat("Error:", e$message, "\n")
    })
  }
}

# perform_DE_compareDays_onCelltypeStrain(allvirus, "outputs/all_virus/DE_compareDays")