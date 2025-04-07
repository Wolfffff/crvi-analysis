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
# both: total read counts (ALTso give  + REF) (same dimensions)
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

# Preallocate only for the current chunk
pvals1 <- data.frame(site = site_ids, matrix(NA, n_chunk, 4))
coefs1 <- data.frame(site = site_ids, matrix(NA, n_chunk, 4))
se1    <- data.frame(site = site_ids, matrix(NA, n_chunk, 4))
lik1   <- data.frame(site = site_ids, logLik_mod1 = NA_real_, logLik_mod2 = NA_real_)
var1   <- data.frame(site = site_ids, var = NA_real_)
#
fit_issues <- data.frame(
  site = site_ids,
  mod1_status = rep(NA_character_, n_chunk),
  mod2_status = rep(NA_character_, n_chunk),
  stringsAsFactors = FALSE
)


# Use names from glmmTMB output summary table
colnames(pvals1)[2:5] <- colnames(coefs1)[2:5] <- colnames(se1)[2:5] <-
  c("(Intercept)", "treatmentControl", "treatment3mM", "sexmale")

# Initialize log for downsampled samples (all logged once per chunk)
# ------------------------------------------------------------------------------
log_chunk <- data.frame(chr_pos = character(),
                        sample = character(),
                        orig_alt = integer(),
                        orig_both = integer(),
                        downsampled_alt = integer(),
                        stringsAsFactors = FALSE)

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
  idx_high_depth <- which(info$tmp_y > 20) # Samples needing downsampling

  # ----------------------------------------------------------------------------
  # Downsample samples with depth > 20 using sampling WITHOUT replacement
  # ----------------------------------------------------------------------------
  if (length(idx_high_depth) > 0) {
    for (j in idx_high_depth) {
      sample_id <- info$sample[j]
      alt_j <- info$tmp_x[j]
      total_j <- info$tmp_y[j]
      ref_j <- total_j - alt_j

      if (total_j >= 20) {
        # sample without replacement using hypergeometric distribution
        down_alt <- rhyper(1, alt_j, ref_j, 20)
        # Update data with downsampled values
        info$tmp_x[j] <- down_alt
        info$tmp_y[j] <- 20

        # Log the downsampling event
        log_chunk <- rbind(log_chunk, data.frame(
          chr_pos = chr_pos,
          sample = sample_id,
          orig_alt = alt_j,
          orig_both = total_j,
          downsampled_alt = down_alt,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # ----------------------------------------------------------------------------
  # Subset to valid samples: must have reads and known sex
  # ----------------------------------------------------------------------------
  tmp_info <- subset(info, tmp_y > 0 & sex != 'unknown')

  # ----------------------------------------------------------------------------
  # Fit beta-binomial regression models
  # Full model includes "condition" (effect of interest)
  # Reduced model omits "condition" for likelihood ratio test
  # ----------------------------------------------------------------------------
  # Full model with treatment
# Fit mod1, initializing status for warnings
mod1_status <- "OK"
mod1_warnings <- character()

#run with warning and error logging
mod1 <- tryCatch({
  withCallingHandlers({
    result <- glmmTMB(
      cbind(tmp_x, tmp_y - tmp_x) ~ treatment + sex + (1 | treatment:plate),
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

# Combine warnings only if no error occurred
if (!startsWith(mod1_status, "Error") && length(mod1_warnings) > 0) {
  mod1_status <<- paste("Warnings:", paste(mod1_warnings, collapse = "; "))
}
# Fit mod2
mod2_status <- "OK"
mod2_warnings <- character()

mod2 <- tryCatch({
  withCallingHandlers({
    result <- glmmTMB(
      cbind(tmp_x, tmp_y - tmp_x) ~ sex + (1 | treatment:plate),
      family = betabinomial(),
      data = tmp_info,
      control = control_list
    )
    result  # this needs to be returned explicitly!
  }, warning = function(w) {
    mod2_warnings <<- c(mod2_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
}, error = function(e) {
  mod2_status <<- paste("Error:", conditionMessage(e))
  NULL
})

# Combine warnings only if no error occurred
if (!startsWith(mod2_status, "Error") && length(mod2_warnings) > 0) {
  mod2_status <<- paste("Warnings:", paste(mod2_warnings, collapse = "; "))
}

  #Log results (loglik)
  tryCatch({ lik1$logLik_mod1[k] <- logLik(mod1) }, error = function(e) {})
  tryCatch({ lik1$logLik_mod2[k] <- logLik(mod2) }, error = function(e) {})
# log results (summary outputs)
  tryCatch({
    fixed <- summary(mod1)$coefficients$cond
    coefs1[k, names(fixed[, 1])] <- fixed[, 1]
    se1[k, names(fixed[, 2])]    <- fixed[, 2]
    pvals1[k, names(fixed[, 4])] <- fixed[, 4]
  }, error = function(e) {})
#log variance associated with the random effect
  tryCatch({
    var1$var[k] <- VarCorr(mod1)$cond[["treatment:plate"]][1, 1]
  }, error = function(e) {})
#log fit issues
  fit_issues$mod1_status[k] <- mod1_status
  fit_issues$mod2_status[k] <- mod2_status


  print(i)
}

# ------------------------------------------------------------------------------
# Write out results for this chunk (filenames include the starting row)
# ------------------------------------------------------------------------------
write.table(pvals1, paste0('/Genomics/argo/users/ed7982/crvi-analysis/longevity_bp_resolution/bbin/output/Apr6_bbin_pvals_', ITER1, '.txt'),
            row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)

write.table(coefs1, paste0('/Genomics/argo/users/ed7982/crvi-analysis/longevity_bp_resolution/bbin/output/Apr6_bbin_coefs_', ITER1, '.txt'),
            row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)

write.table(se1, paste0('/Genomics/argo/users/ed7982/crvi-analysis/longevity_bp_resolution/bbin/output/Apr6_bbin_se_', ITER1, '.txt'),
            row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)

write.table(lik1, paste0('/Genomics/argo/users/ed7982/crvi-analysis/longevity_bp_resolution/bbin/output/Apr6_bbin_loglik_', ITER1, '.txt'),
            row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)

write.table(var1, paste0('/Genomics/argo/users/ed7982/crvi-analysis/longevity_bp_resolution/bbin/output/Apr6_bbin_var_', ITER1, '.txt'),
            row.names = FALSE, sep = '\t', col.names = TRUE, quote = FALSE)

write.table(fit_issues,
  file = paste0('/Genomics/argo/users/ed7982/crvi-analysis/longevity_bp_resolution/bbin/output/Apr6_bbin_fit_issues_', ITER1, '.txt'),
  sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)


if (nrow(log_chunk) > 0) {
  write.table(log_chunk,
              file = paste0('/Genomics/argo/users/ed7982/crvi-analysis/longevity_bp_resolution/bbin/output/Apr6_bbin_downsampled_sites', ITER1, '.tsv'),
              sep = '\t',
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE)
}
