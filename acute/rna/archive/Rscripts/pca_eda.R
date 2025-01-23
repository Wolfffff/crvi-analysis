rnaseq_data <- import("cache/rnaseq_all_2024-02-07.rds")
covariates  <- import("cache/covariates_2024-02-07.rds")

source(here::here("rna/Rscripts/functions.R"))


pca_and_plot <- function(data, covariates, column_name, output_file = "pca_plot.jpg") {
  # Perform PCA
  pca_res <- pca(data)
  
  # Create a data frame for plotting
  pca_df <- data.frame(t(pca_res$pc))
  pca_df[['plate']] <- covariates[['plate']]
  pca_df[['treatment']] <- covariates[['treatment']]
  
  # Plot + save

  ggplot(pca_df, aes(x = PC1, y = PC2, shape=treatment, color=plate)) +
    geom_point() +
    labs(color = column_name) +
    ggtitle("PCA Plot") +
    theme_minimal()

  ggsave(output_file, width = 10, height = 10, dpi = 300)
}

pca_and_plot(rnaseq_data$body$cpm, rnaseq_data$body$covariates, 'plate', 'pca_plot_raw.jpg')
pca_and_plot(rnaseq_data$body$batch_residuals$cpm, rnaseq_data$body$covariates, 'plate', 'pca_batch_resid.jpg')
