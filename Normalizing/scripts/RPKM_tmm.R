#
# Calculating TMM-normalized RPKM-values
#
source("scripts/Functions.R")
suppressPackageStartupMessages(library(edgeR))
save.image(paste0("images/", snakemake@rule, ".RData"))

for(i in 1:length(snakemake@input)){
  load(snakemake@input[[i]])

  RPKM = calc_rpkm(txi$counts, txi$length)
  RPKM = round(RPKM)
  
  y = prep.DGEList(txi, RPKM)
  txi$counts = calc_tmm(y)
  
  save(txi, file = snakemake@output[[i]])
}