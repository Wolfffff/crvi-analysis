source(here::here("rna/Rscripts/load_install_packages.R"))
source(here::here("rna/Rscripts/functions.R"))

countdata = list(head = read_tsv(here::here("rna/head/results/counts/all.symbol.tsv")),
                 body = read_tsv(here::here("rna/body/results/counts/all.symbol.tsv")))

covariates = rbind(read_tsv(here::here("rna/body/config/samples.tsv")),
                   read_tsv(here::here("rna/head/config/samples.tsv")))

filterGenes <- function(countdata.norm){
    print(paste("Total genes:", nrow(countdata.norm$genes)))
    # Removing genes where max expression is less than 1 cpm
    max_cpm = apply(cpm(countdata.norm), 1, max)
    drop <- which(max_cpm < 1)
    print(length(drop))
    countdata.norm <- countdata.norm[-drop, ]

    # Removing genes with average expression less than 1 cpm
    mean_cpm = apply(cpm(countdata.norm), 1, mean)
    drop <- which(mean_cpm < 1)
    print(length(drop))
    countdata.norm <- countdata.norm[-drop, ]

    # Removing genes with mean expression less than 10 counts
    mean_count = apply(countdata.norm$counts, 1, mean)
    drop <- which(mean_count < 3)
    print(length(drop))
    countdata.norm <- countdata.norm[-drop, ]

    # Removing genes with expression less than 1 cpm is more than 20% of the samples
    prop_non_expressed = apply(cpm(countdata.norm), 1, \(x) sum(x < 1)) / nrow(countdata.norm$samples)
    drop <- which(prop_non_expressed > 0.2)
    print(length(drop))
    countdata.norm <- countdata.norm[-drop, ]

    print(paste("Filtered genes:", nrow(countdata.norm$genes)))
    countdata.norm
}

print("Filtering genes")
y <- map(countdata, DGEList)
y <- map(y, calcNormFactors)
y.filtered <- map(y, filterGenes)


#################
# Set model matrices
#################

correctModelMatrix = function(X){
    QR <- qr(crossprod(X))                 # Get the QR decomposition
    vars <- QR$pivot[seq_len(QR$rank)]     # Variable numbers
    X = X[,vars]
    X
}

setModelMatrices = function(current_tissue, y, covariates){
    cov = filter(covariates, tissue == current_tissue) |> 
          filter(sample_name %in% rownames(y[[current_tissue]]$samples)) |>
          arrange(match(sample_name, rownames(y[[current_tissue]]$samples)))
    mod1 <- model.matrix(~1 + plate + treatment, cov) 
    mod0 <- model.matrix(~1 + plate, cov) 
    design <- model.matrix(~as.factor(treatment), cov) 
    colnames(design)[1:2] <- c("Intercept","crvi")
    colnames(mod1)[1:2] <- c("Intercept","crvi")
    rownames(y[[current_tissue]]$counts) <- y[[current_tissue]]$genes[,1]
    list(counts = y[[current_tissue]], 
         tissue = current_tissue,
         covariates = cov, 
         mod1 = correctModelMatrix(mod1),
         mod0 = correctModelMatrix(mod0), 
         design = design)
}

print("Setting model matrices")
rnaseq_data = map(c(head = "head", body = "body"), setModelMatrices, y.filtered, covariates)

setCPMVoom <- function(data){
    data$cpm  <-  cpm(data$counts, log=TRUE)
    data$voom <- voom(data$counts, data$mod1)
    data$l2c  <- log2(data$counts$counts + 1)
    data
}

setSVAnum = function(data){
    counts_list = list(cpm = data$cpm, 
                       voom = data$voom$E, 
                       l2c = data$l2c)
    # Only if n
    data$n.sva <- future_map(counts_list, num.sv, data$mod1, method="leek")
    data
}
print("Setting CPM voom")
rnaseq_data = map(rnaseq_data, setCPMVoom)

print("Setting SVA")
rnaseq_data = future_map(rnaseq_data, setSVAnum)
rnaseq_data |> map("n.sva")

setSVA = function(data){
     counts_list = list(cpm = data$cpm, 
                        voom = data$voom$E, 
                        l2c = data$l2c)
     
     data$sva <- future_map2(counts_list[data$n.sva!=0], data$n.sva[data$n.sva!=0],
                             \(x, y) sva(x, data$mod1, data$mod0, n.sv=y, method = "two-step"))
     data
}
rnaseq_data = map(rnaseq_data, setSVA)

getBatchResiduals <- function(x, label, tissue, design, mod0, sva = NULL, treatment){
    no.batch <- removeBatchEffect(x, 
                                  design = design, 
                                  covariates = cbind(mod0, sva))
    pca_no.batch = pca(no.batch)
    pdf(paste0('tmp/PCA-batch-corrected-', tissue, '-', label, '.pdf'), width = 1080, height = 1080)
        print(pca_plot(pca_no.batch, c("C", "CrVI")[treatment]))
    dev.off()
    rownames(no.batch) = data$counts$genes[,1]
    no.batch
}

makeResiduals <- function(data){
    covariates <- data$covariates
    counts_list = list(cpm = data$cpm, 
                       voom = data$voom$E, 
                       l2c = data$l2c)
    data$batch_residuals <- future_map2(counts_list, names(counts_list), 
                            \(x, y) getBatchResiduals(x, y, data$tissue, 
                                                      data$design, data$mod0,
                                                      treatment = as.numeric(as.factor(pull(covariates, treatment)))))  
    data
}

# Futures max memory
options(future.globals.maxSize = 1e9 * 1024)
future::plan(list(future::tweak(future::multisession, workers = 2), 
                  future::tweak(future::multisession, workers = 3)))

rnaseq_data = future_map(rnaseq_data, makeResiduals)

setMouthwash = function(data){
     counts_list = list(cpm = data$cpm, 
                        voom = data$voom$E, 
                        l2c = data$l2c)
     
     data$mwash <- future_map2(counts_list[data$n.sva!=0], data$n.sva[data$n.sva!=0], 
                             \(x, y) mouthwash(Y = t(x), X = data$mod1, k = y, cov_of_interest = 2,
                                               include_intercept = FALSE))
     data
}
rnaseq_data = future_map(rnaseq_data, setMouthwash)

makeResidualsMwash <- function(data){
    covariates <- data$covariates
    counts_list = list(cpm = data$cpm, 
                       voom = data$voom$E, 
                       l2c = data$l2c)
    data$mwash_residuals <- future_map2(counts_list[data$n.sva!=0], names(counts_list)[data$n.sva!=0], 
                            \(x, y) getBatchResiduals(x, paste0(y, "-mwash"), data$tissue, 
                                                      data$design, data$mod0, sva = data$mwash[[y]]$Zhat,
                                                      treatment = as.numeric(as.factor(pull(covariates, treatment)))))  
    data
}
rnaseq_data = future_map(rnaseq_data, makeResidualsMwash)

export(rnaseq_data, affix_date("cache/rnaseq_all.rds"))
export(covariates, affix_date("cache/covariates.rds"))

# rnaseq_data <- import("rna/cache/rnaseq_all_2023-12-04.rds")
# covariates  <- import("rna/cache/covariates_2023-12-04.rds")
