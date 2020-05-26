#
# Calculating cpm values with loess normalization (voom)
#
source("scripts/Functions.R")
save.image(paste0("images/", snakemake@rule, ".RData"))
suppressPackageStartupMessages(library(limma))

for(i in 1:length(snakemake@input)){
  load(snakemake@input[[i]])
  counts = txi$counts
  cl = limma::normalizeCyclicLoess(counts)
  cl[which(cl < 0)] = 0
  txi$counts = cl
  save(txi, file = snakemake@output[[i]])
}