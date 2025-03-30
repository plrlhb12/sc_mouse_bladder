# still DE but compare between ctrl and virus samples on each day for each cell type; 
perform_DE_compareCtrlVirus_onTypeDay <- function(all.obj) {
  # set variables
  whole.obj <- all.obj
  DefaultAssay(whole.obj) <- "RNA"
  cell.types <- levels(factor(whole.obj$cell.type))
  days <- levels(factor(whole.obj$days))
  sars.status <- levels(factor(whole.obj$groups))
  
  # start day loop
  Idents(whole.obj) <- "days"
  for (day in days){
    day.obj <- subset(whole.obj, idents = day)
    parent_dir <- paste0("outputs/", day, "_allStrainTogether")# here may need to customize
    dir.create(parent_dir) # here may need to customize
    
    Idents(day.obj) <- "cell.type"
    # start cell type loop to compare sars ctrl VS infected
    for (type in cell.types){
      print(paste0("looping ", type, "for ", day))
      # add error handle
      tryCatch({
        sub.obj <- subset(day.obj, idents = type)
        run.name <- type
        final.out.dir <- paste0(parent_dir, "/", run.name)
        dir.create(final.out.dir)
        
        # set compare target of interest, here is to compare groups status
        compare_list <- sars.status
        print(paste0("compare: ", compare_list))
        num_objects <- length(compare_list)
        
        Idents(sub.obj) <- "groups"
        for (i in 1:(num_objects - 1)) {
          for (j in (i + 1):num_objects) {
            if (i != j) {
              compare <- paste0(compare_list[[i]], "vs", compare_list[[j]])
              print(compare)
              ident1 = compare_list[[i]]
              ident2 = compare_list[[j]]
              
              # Use tryCatch to handle the error
              tryCatch({
                markers <- FindMarkers(sub.obj, ident.1 = ident1, ident.2 = ident2)
                write.csv(markers, file = paste0(final.out.dir, "/", run.name, "_", compare, "_DE_genes.csv"))
                
                top50 <- markers %>%
                  slice_max(n = 50, order_by = avg_log2FC)
                write.csv(top50, file = paste0(final.out.dir, "/", run.name, "_", compare, "_top50_genes.csv"))
                
                bottom50 <- markers %>%
                  slice_min(n = 50, order_by = avg_log2FC)
                write.csv(bottom50, file = paste0(final.out.dir, "/", run.name, "_", compare, "_bottom50_genes.csv"))
                
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
}

# source perform_DE_compareCtrlVirus_onTypeDay.R
# perform_DE_compareCtrlVirus_onTypeDay(all.obj)