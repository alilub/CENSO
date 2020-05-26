#
# Calculating TMM-normalized TPM values
#
source("scripts/Functions.R")
suppressPackageStartupMessages(library(edgeR))
save.image(paste0("images/", snakemake@rule, ".RData"))

for(i in 1:length(snakemake@input)){
  load(snakemake@input[[i]])
  
  tpm = calc_tpm(txi$counts, txi$length) 
  norm.factors = calcNormFactors(tpm)
  lib.size = colSums(tpm)
  
  tmm_scaleFactors <- lib.size * norm.factors
  tmm_normFactors <- tmm_scaleFactors/exp(mean(log(tmm_scaleFactors)))
  
  counts_tmm <- tpm %*% diag(1/tmm_normFactors)
  colnames(counts_tmm) <- colnames(tpm)
  
  counts_tmm <- round(counts_tmm)
  sf = 1/tmm_normFactors
  names(sf) = colnames(txi$counts)
  save(sf, file = snakemake@output[[i+length(snakemake@input)]])
  
  txi$counts = counts_tmm
  save(txi, file = snakemake@output[[i]])
}