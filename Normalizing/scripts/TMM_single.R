#
# Calculating TMM-normalized counts
#
source("scripts/Functions.R")
suppressPackageStartupMessages(library(edgeR))

save.image(paste0("images/", snakemake@rule, ".RData"))

for(i in 1:length(snakemake@input)){
  load(snakemake@input[[i]])
  y = prep.DGEList(txi = txi)
  txi$counts = calc_tmm(y)
  save(txi, file = snakemake@output[[i]])
}