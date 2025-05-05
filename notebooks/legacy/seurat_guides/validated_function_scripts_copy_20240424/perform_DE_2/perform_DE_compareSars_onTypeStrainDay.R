# still DE but compare between POS and NEG for each cell type of each strain on each day
perform_DE_compareSars_onTypeStrainDay <- function(allvirus) {
  # set variables
  whole.obj <- allvirus
  DefaultAssay(whole.obj) <- "RNA"
  cell.types <- levels(factor(whole.obj$cell.type))
  strains <- levels(factor(whole.obj$strains))
  days <- levels(factor(whole.obj$days))
  sars.status <- levels(factor(whole.obj$sars.pos))
  
  # start day loop
  Idents(whole.obj) <- "days"
  for (day in days){
    day.obj <- subset(whole.obj, idents = day)
    parent_dir <- paste0("outputs/", day)
    dir.create(parent_dir)
    
    # error handler since some samples don't have D10
    tryCatch({
      # start strain loop
      Idents(day.obj) <- "strains"
      for (strain in strains){
        strain.day.obj <- subset(day.obj, idents = strain)
        dir.create(paste0(parent_dir, "/", strain))
        
        # start cell type loop to compare sars POS VS NEG
        print(paste0("start looping cell type for ", strain))
        Idents(strain.day.obj) <- "cell.type"
        for (type in cell.types){
          print(paste0("looping ", type, " for ", strain))
          # add error handle
          tryCatch({
            sub.obj <- subset(strain.day.obj, idents = type)
            run.name <- type
            final.out.dir <- paste0(parent_dir, "/", strain, "/", run.name)
            dir.create(final.out.dir)
            
            # set compare target of interest, here is to compare sars status
            compare_list <- sars.status
            print(paste0("compare: ", compare_list))
            num_objects <- length(compare_list)
            
            Idents(sub.obj) <- "sars.pos"
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
      }
      
      # ... (rest of the code for writing CSV files)
    }, error = function(e) {
      # Print a message or take other actions upon error
      cat("Error:", e$message, "\n")
    })
  }
}

# perform_DE_compareSars_onTypeStrainDay(allvirus)