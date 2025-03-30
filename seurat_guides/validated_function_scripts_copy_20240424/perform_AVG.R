#############################################################
# using AverageExpression() function to calculate average expression for each cluster
#############################################################
# by default, all the applicable assays including RNA, refAssay, integrated, and etc., will perform the function at the same time
# the returned object e.g., avg.exp is a list of matrix instead of a single matrix
# according to the active idents, and return a matrix for each assay
# if not define "slot = ", the default slot is "data", i.e., log transformed and normalized data
# The ideal assay for calculating AVG gene expression should be the one from "RNA" assay; do not use "integrated" for expression analysis; "SCT" depends.

Idents(d3) <- d3$seurat_clusters
avg.exp <- AverageExpression(d3) # calculaiton AVG on the slot of "data" for all the assays and export a list of matrix
for (assay_id in names(avg.exp)) {
  avg.exp.df <- as.data.frame(avg.exp[[assay_id]])
  write.csv(avg.exp.df, file = paste0("d3", assay_id, "_avg_exp.csv"), row.names = TRUE)
}

avg.exp <- AverageExpression(d3,slot = "counts") # if want to explore using counts slot
for (assay_id in names(avg.exp)) {
  avg.exp.df <- as.data.frame(avg.exp[[assay_id]])
  write.csv(avg.exp.df, file = paste0("d3", assay_id, "_avg_counts.csv"), row.names = TRUE)
}

for (i in batches){
  sc <- seurat_list[[i]]
  Idents(sc) <- sc$predicted.celltype
  avg.exp <- as.data.frame(AverageExpression(sc, assays = "refAssay")) # only peform AVG on 1 assay
  write.csv(avg.exp, file = paste0(run.name, "/", i, "_avg_exp_by_type.csv"))
  }

# most important things: define active assay using "RNA" inside the funciton
DefaultyAssay(obj) <- "RNA" # this doesn't help; need to define assays = "RNA" inside the function if only want to perform calculation on the "RNA" assay
Idents(sub.obj) <- sub.obj$predicted.celltype # the active.idents actually doesn't affect the calculation
avg.exp <- as.data.frame(AverageExpression(sub.obj, assays = "RNA", group.by = c("days"))) # groups by 1 conditon
write.csv(avg.exp, file = paste0(run.name, "/gene_avg_exp_by_days.csv"))
avg.exp <- as.data.frame(AverageExpression(sub.obj, assays = "RNA", group.by = c("strains")))
write.csv(avg.exp, file = paste0(run.name, "/gene_avg_exp_by_strains.csv"))
avg.exp <- as.data.frame(AverageExpression(sub.obj, assays = "RNA", group.by = c("strains", "days"))) # groups by 2 conditon
write.csv(avg.exp, file = paste0(run.name, "/gene_avg_exp_by_strains_days.csv"))
avg.exp <- as.data.frame(AverageExpression(sub.obj, assays = "RNA", group.by = c("sars.pos")))
write.csv(avg.exp, file = paste0(run.name, "/gene_avg_exp_by_sars.csv"))
avg.exp <- as.data.frame(AverageExpression(sub.obj, assays = "RNA", group.by = c("is.Omicron")))
write.csv(avg.exp, file = paste0(run.name, "/gene_avg_exp_OmicronVSothers.csv"))