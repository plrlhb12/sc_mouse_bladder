# here is to use rename function in R to rename all files by replacing some strings
directory_path <- "../20230823-individual/individuals-2"
# make changes on pattern = "newpattern"
files <- list.files(directory_path, full.names = TRUE, pattern = "\\.rds$")

# generate a list of 8 sublists. Each sublist has two elements. 
# The name of the first element is "from", the second is "to"
# 1.change D3 and D7 to D03 and D07, respectively
# 2.change mu to Mu
# 3.change delta to Delta; Delta-Plus to Delta_Plus
# 4.rename Omicron to reflect the clade
# 5.remove first 11 characters "2023-08-28/29"

patterns <- list(
  list(from = "^...........", to = ""),
  list(from = "D3", to = "D03"),
  list(from = "D7", to = "D07"),
  list(from = "mu", to = "Mu"),
  list(from = "delta", to = "Delta"),
  list(from = "4J-Omicron", to = "4J-Omicron-B.1.1.529"),
  list(from = "4K-Omicron", to = "4K-Omicron-B.1.529"),
  list(from = "4M-Omicron", to = "4M-Omicron-BA.2")
)

# Function to apply renaming patterns
rename_files <- function(files, patterns) {
  #loop each filename
  for (file in files) {
    new_name <- basename(file)
    # loop each pattern for each filename
    for (pattern in patterns) {
      from_pattern <- pattern$from # get the value of the element of "from"
      to_pattern <- pattern$to
      new_name <- gsub(from_pattern, to_pattern, new_name)
    }
    
    new_path <- file.path(dirname(file), new_name)
    file.rename(file, new_path)
    print(paste("Renamed:", basename(file), "to", new_name))
  }
}

rename_files(files, patterns)
print(files)