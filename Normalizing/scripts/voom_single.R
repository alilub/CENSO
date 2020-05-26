#
# Calculating cpm values with loess normalization (voom)
#
source("scripts/Functions.R")
save.image(paste0("images/", snakemake@rule, ".RData"))
suppressPackageStartupMessages(library(limma))

for(i in 1:length(snakemake@input)){
  load(snakemake@input[[i]])
  counts = txi$counts
  Rj = colSums(counts)
  for(j in 1:length(Rj)){
    counts[,j] = ((counts[,j] + 0.5)/(Rj[j] + 1.0))*10^6
  }
  cl = affy::normalize.loess(counts, log.it = F, family.loess = "gaussian")
  cl[which(cl < 0)] = 0
  txi$counts = cl
  save(txi, file = snakemake@output[[i]])
}