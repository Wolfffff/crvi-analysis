library(rio)
library(tidyverse)

# rna_seq_file = "cache/rnaseq_all_2023-12-04.rds"
rna_seq_file = snakemake@input[["rna_seq"]]
all_data = import(rna_seq_file)

sample_list_file = snakemake@input[["sample_list"]]

tissue = snakemake@wildcards[["tissue"]]
data = all_data[[tissue]]

samples = read.table(sample_list_file, header = FALSE, sep = "\t")[,1]

colnames(data$mwash$l2c$Zhat) = paste0("Zhat", 1:dim(data$mwash$l2c$Zhat)[2])

# Load map and get dna from 
if (tissue == "head") {
  mapping_df  <-  read.csv("/Genomics/argo/users/swwolf/git/crvi-acute/acute_metdata.csv") 
  # map sample_name rna_head and get dna from row
  matched_df <- data$covariates |>
    left_join(mapping_df, by = c("sample_name" = "rna_head"))
  data$covariates$sample_name = matched_df$dna
}


input_df = cbind(data$covariates, 
                 data$mwash$l2c$Zhat) |>
               filter(sample_name %in% samples)


input_df = input_df[match(samples, input_df$sample_name),]

stopifnot(all(input_df$sample_name == samples))
# Go back to original names
samples = all_data[[tissue]]$covariates$sample_name

phenos = data$l2c[,samples]

write.table(rownames(phenos), snakemake@output[["gene_list"]], 
          row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
export(input_df, snakemake@output[["covariates"]])
export(phenos, snakemake@output[["phenotypes"]])

