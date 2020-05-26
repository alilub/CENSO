#
# Caclulating DESeq-normalized TPM values
#

#
# DESeq normalzing method
#
norm_factors <- function(mat) {
  nz <- apply(mat, 1, function(row) !any(round(row) == 0))
  mat_nz <- mat[nz, , drop = FALSE]
  p <- ncol(mat)
  geo_means <- exp(apply(mat_nz, 1, function(row) mean(log(row))))
  s <- sweep(mat_nz, 1, geo_means, `/`)
  sf <- apply(s, 2, median)
  scaling <- exp( (-1 / p) * sum(log(sf)))
  sf * scaling
}

source("scripts/Functions.R")
save.image(paste0("images/", snakemake@rule, ".RData"))
suppressPackageStartupMessages(library(DESeq2))
for(i in 1:length(snakemake@input)){
    load(snakemake@input[[i]])
    #
    # Normalization
    #
    tpm = calc_tpm(txi$counts, txi$length)
    tpm_sf <- norm_factors(tpm)
    tpm_norm <- as.data.frame(t(t(tpm) / tpm_sf))
    
    txi$counts = as.matrix(tpm_norm)
    save(txi, file = snakemake@output[[i]])
    sf = 1/tpm_sf
    save(sf, file = snakemake@output[[i+length(snakemake@input)]])
}

