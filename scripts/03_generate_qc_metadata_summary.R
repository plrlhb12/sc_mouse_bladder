library(Seurat)
library(dplyr)
library(readr)
library(here)

# Load your unfiltered and filtered Seurat data
seurat_unfiltered <- readRDS(here("results", "seurat_list_unfiltered.rds"))
merged_filtered <- readRDS(here("results", "merged_filtered_percentile_filtered.rds"))

# Split merged_filtered back into per-sample objects (if needed)
filtered_list <- SplitObject(merged_filtered, split.by = "sample")

# Initialize storage
all_meta <- list()
summary_stats <- list()

# Compare each sample
for (sample in names(seurat_unfiltered)) {
  unfiltered <- seurat_unfiltered[[sample]]
  filtered <- filtered_list[[sample]]
  
  # Add percent.mt if not already there
  if (!"percent.mt" %in% colnames(unfiltered@meta.data)) {
    unfiltered[["percent.mt"]] <- PercentageFeatureSet(unfiltered, pattern = "^mt-")
  }
  if (!"percent.mt" %in% colnames(filtered@meta.data)) {
    filtered[["percent.mt"]] <- PercentageFeatureSet(filtered, pattern = "^mt-")
  }
  
  # Metadata
  meta_pre <- unfiltered@meta.data %>%
    mutate(cell = colnames(unfiltered), sample = sample, status = "Pre-filter")
  
  meta_post <- filtered@meta.data %>%
    mutate(cell = colnames(filtered), sample = sample, status = "Post-filter")
  
  all_meta[[sample]] <- bind_rows(meta_pre, meta_post)
  
  # Summary stats
  summary_stats[[sample]] <- tibble(
    sample = sample,
    total_cells = ncol(unfiltered),
    kept_cells = ncol(filtered),
    removed_cells = ncol(unfiltered) - ncol(filtered),
    percent_retained = round(100 * ncol(filtered) / ncol(unfiltered), 1)
  )
}

# Combine all and save
qc_metadata <- bind_rows(all_meta)
qc_summary <- bind_rows(summary_stats)

write_tsv(qc_metadata, here("results", "qc_metadata_before_after.tsv"))
write_tsv(qc_summary, here("results", "qc_filtering_summary.tsv"))

message("✅ QC metadata and summary written to 'results/'")

