#
# Calculating RPKM values
#
source("scripts/Functions.R")
suppressPackageStartupMessages(library(limma))
save.image(paste0("images/", snakemake@rule, ".RData"))

for(i in 1:length(snakemake@input)){
  print(snakemake@input[[i]])
  load(snakemake@input[[i]])
  counts = calc_rpkm(txi$counts, txi$length)
  cl = limma::normalizeCyclicLoess(counts)
  cl[which(cl < 0)] = 0
  txi$counts = cl
  save(txi, file = snakemake@output[[i]])
}