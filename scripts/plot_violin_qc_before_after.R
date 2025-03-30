library(ggplot2)
library(dplyr)
library(readr)
library(here)

# Define the plotting function
plot_violin_comparison <- function(df, metric, y_label, filename) {
  # Create "Spacer" rows between Pre-filter and Post-filter
  spacer_df <- df %>%
    distinct(sample) %>%
    mutate(status = "Spacer", !!metric := NA_real_)
  
  # Combine spacer and original data
  df_spaced <- bind_rows(df, spacer_df) %>%
    mutate(status = factor(status, levels = c("Pre-filter", "Spacer", "Post-filter")))
  
  # Create violin plot
  p <- ggplot(df_spaced, aes(x = sample, y = .data[[metric]], fill = status)) +
    geom_violin(data = filter(df_spaced, status != "Spacer"),
                alpha = 0.5,
                position = position_dodge(width = 0.8)) +
    geom_jitter(data = filter(df_spaced, status != "Spacer"),
                aes(color = status),
                alpha = 0.2,
                size = 0.2,
                position = position_jitterdodge(jitter.width = 0.2)) +
    coord_flip() +
    labs(x = NULL, y = y_label, title = paste("Before vs After Filtering:", y_label)) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.position = "top"
    ) +
    scale_fill_manual(values = c("Pre-filter" = "skyblue", "Post-filter" = "forestgreen")) +
    scale_color_manual(values = c("Pre-filter" = "skyblue", "Post-filter" = "forestgreen"))
  
  # Save plot
  ggsave(filename = here("results", filename), plot = p, width = 9, height = 5)
}

# Create plots for each QC metric
plot_violin_comparison(meta_all, "nFeature_RNA", "# of Detected Genes", "violin_nFeature_before_after.pdf")
plot_violin_comparison(meta_all, "percent.mt", "Mitochondrial Percentage", "violin_percent_mt_before_after.pdf")
plot_violin_comparison(meta_all, "nCount_RNA", "UMI Count per Cell", "violin_nCount_before_after.pdf")
