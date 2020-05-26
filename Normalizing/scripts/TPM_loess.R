#
# Calculating TMM-normalized counts
#
source("scripts/Functions.R")
save.image(paste0("images/", snakemake@rule, ".RData"))

for(i in 1:length(snakemake@input)){
  load(snakemake@input[[i]])
  tpm = calc_tpm(txi$counts, txi$length)
  cl = limma::normalizeCyclicLoess(tpm)
  cl[which(cl < 0)] = 0
  txi$counts = cl
  save(txi, file = snakemake@output[[i]])
}