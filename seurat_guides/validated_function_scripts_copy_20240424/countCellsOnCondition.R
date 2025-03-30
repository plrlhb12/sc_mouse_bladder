# calculate the frequency of the categories in a column, e.g., calculate cell number per type
table(obj@meta.data$days)

# if want to save the std into a text file
sink("allvirus_sars_cellType_day.txt")
frequency_table <- table(meta$cell.type, meta$days, meta$sars.pos)
print(frequency_table)
sink()

# if want to save the std out to an execel file and covert a long form to a wide form
a.table <- table(sub.obj$predicted.celltype, sub.obj$days)
df <- data.frame(a.table)
df <- pivot_wider(df, id_cols = Var1, names_from = Var2, values_from = Freq)
write.csv(df, file = paste0(run.name, "/", run.name, "_cell_numbers_per_type.csv"), row.names = FALSE)

# if want to perform calculation on multiple objects
# generate cell numbers by type according to conditions
conditions <- c("sars.pos", "strains", "days", "replicates", "is.Omicron")
for (condition in conditions){
  a.table <- table(sub.obj$predicted.celltype, sub.obj@meta.data[[condition]])
  df <- data.frame(a.table)
  df <- pivot_wider(df, id_cols = Var1, names_from = Var2, values_from = Freq)
  write.csv(df, file = paste0("all_virus/all_virus_s_", condition, ".csv"), row.names = FALSE)
