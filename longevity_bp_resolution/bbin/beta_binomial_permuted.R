#!/usr/bin/env Rscript

#import libraries, seed, parameters for GLMM
library(data.table)
library(glmmTMB)
set.seed(123)
control_list <- glmmTMBControl(optCtrl = list(iter.max = 2e3, eval.max = 2e3))

#read in start from command line and get range (either do 5000 sites or finish, whichever comes first)
args <- commandArgs(trailingOnly = TRUE)
ITER1 <- as.integer(args[1])
ITER2 <- min(ITER1 + 4999, 1119387)


# ------------------------------------------------------------------------------
# Load input data:
# alt: ALT allele counts (sites x samples)
# both: total read counts (same dimensions)
# info: sample metadata (samples x covariates)
# ------------------------------------------------------------------------------
alt = fread('/Genomics/argo/users/ed7982/crvi-analysis/longevity_bp_resolution/bbin/input/allele_table_alt_bbin_with_header_goodplates_fixedmissing.txt')
both = fread('/Genomics/argo/users/ed7982/crvi-analysis/longevity_bp_resolution/bbin/input/allele_table_both_bbin_with_header_goodplates.txt')
info = read.delim('/Genomics/argo/users/ed7982/crvi-analysis/longevity_bp_resolution/bbin/input/sample_names_info.txt')

# ------------------------------------------------------------------------------
# Sanity checks to ensure sample IDs match across files
# ------------------------------------------------------------------------------
stopifnot(identical(names(alt)[-1], names(both)[-1]))
stopifnot(identical(as.character(info$sample), names(both)[-1]))

# ------------------------------------------------------------------------------
# Initialize storage for regression output
# ------------------------------------------------------------------------------
site_ids <- alt[ITER1:ITER2, 1][[1]]
n_chunk <- length(site_ids)

# For each chunk, we want 10 sites for p-values for each parameter, and a flag for fit issues, but we don't need anything else
ctrl_pvals <- data.frame(site = site_ids, matrix(NA, n_chunk, 10))
chromium_pvals <- data.frame(site = site_ids, matrix(NA, n_chunk, 10))

#
fit_issues <- data.frame(
  site = site_ids,
  matrix(NA, n_chunk, 10),
  stringsAsFactors = FALSE
)


# name colums appropriately
colnames(ctrl_pvals)[2:11] <- c("perm1", "perm2", "perm3", "perm4", "perm5",
                             "perm6", "perm7", "perm8", "perm9", "perm10") 
colnames(fit_issues)[2:11] <- colnames(chromium_pvals)[2:11] <- colnames(ctrl_pvals)[2:11]

# read in downsampling performed by first script for whole chunk
# ------------------------------------------------------------------------------
downsampled_whole_chunk = read.delim(paste0('/Genomics/argo/users/ed7982/crvi-analysis/longevity_bp_resolution/bbin/output/Apr6_bbin_downsampled_sites', ITER1, '.tsv'))

#IMPORTANT: polarize with respect to t0
info$treatment = factor(info$treatment, levels = c("T0", "Control", "3mM"))
# ------------------------------------------------------------------------------
# Main loop: process one site at a time
# ------------------------------------------------------------------------------
for (i in ITER1:ITER2) {
  #need an index to write to the right rows, k is relative position in the chunk but i is absolute index
  k <- i - ITER1 + 1

  info$tmp_x <- as.vector(t(alt[i, -1]))    # ALT counts at site i
  info$tmp_y <- as.vector(t(both[i, -1]))   # Total read counts at site i

  chr_pos <- alt[i, 1][[1]]                 # e.g., "chr2:107391"

  #get downsampling at this specific site
  downsampled_site = downsampled_whole_chunk[downsampled_whole_chunk$chr_pos == chr_pos, ]
  
  #for each indiv downampled at this site, replace alt and ref with the downsampled values
  for (m in 1:nrow(downsampled_site)) {
    new_alt <- downsampled_site$downsampled_alt[m]
    new_both <- 20

    #for matching sample in row m, replace alt$tmp_x with new_alt
    info$tmp_x[info$sample == downsampled_site$sample[m]] <- new_alt
    info$tmp_y[info$sample == downsampled_site$sample[m]] <- new_both

  }
  rm(downsampled_site)
  

  # ----------------------------------------------------------------------------
  # Subset to valid samples: must have reads and known sex
  # ----------------------------------------------------------------------------
  tmp_info <- subset(info, tmp_y > 0 & sex != 'unknown')

  #for each permutation:
  for (permutation in 1:10){
    # Randomly permute the treatment labels
    tmp_info$treatment_permuted <- sample(tmp_info$treatment)

    #initialize status for warnings
    mod1_status <- "OK"
    mod1_warnings <- character()
   
   #run the model with warning and error logging
    mod1 <- tryCatch({
      withCallingHandlers({
        result <- glmmTMB(
          cbind(tmp_x, tmp_y - tmp_x) ~ treatment_permuted + sex + (1 | treatment:plate),
          family = betabinomial(),
          data = tmp_info,
          control = control_list
        )
        result  # must return result explicitly
      }, warning = function(w) {
        mod1_warnings <<- c(mod1_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    }, error = function(e) {
      mod1_status <<- paste("Error:", conditionMessage(e))
      NULL
    })
#write the output of this model (for this site, for this permutation iteration) for both the chromium pval, the control pval, and the fit issues
tryCatch({
  smry <- summary(mod1)
  ctrl_pvals[k, permutation + 1] <- smry$coefficients$cond["treatment_permutedControl", "Pr(>|z|)"]
  chromium_pvals[k, permutation + 1] <- smry$coefficients$cond["treatment_permuted3mM", "Pr(>|z|)"]
}, error = function(e) {
  if (fit_issues[k, permutation + 1] == "OK") {
    fit_issues[k, permutation + 1] <- paste("SummaryError:", conditionMessage(e))
  }
})


# Combine warnings only if no error occurred
if (!startsWith(mod1_status, "Error") && length(mod1_warnings) > 0) {
  mod1_status <<- paste("Warnings:", paste(mod1_warnings, collapse = "; "))
}

#if mod status is ok, write ok, if not write it
if (mod1_status == "OK") {
  fit_issues[k, permutation + 1] <- "OK"
} else {
  fit_issues[k, permutation + 1] <- mod1_status
}

  }

#print individual we've finished
  print(i)
}

#write outputs
write.table(ctrl_pvals, paste0('/Genomics/argo/users/ed7982/crvi-analysis/longevity_bp_resolution/bbin/perm_output/Apr7_bbin_ctrl_pvals_', ITER1, '.txt'),
            row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)

write.table(chromium_pvals, paste0('/Genomics/argo/users/ed7982/crvi-analysis/longevity_bp_resolution/bbin/perm_output/Apr7_bbin_chrom_pvals_', ITER1, '.txt'),
            row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)

write.table(fit_issues,
  file = paste0('/Genomics/argo/users/ed7982/crvi-analysis/longevity_bp_resolution/bbin/perm_output/Apr7_bbin_perm_fit_issues_', ITER1, '.txt'),
  sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
