#
# Calculating RPKM values
#
source("scripts/Functions.R")
save.image(paste0("images/", snakemake@rule, ".RData"))

for(i in 1:length(snakemake@input)){
  print(snakemake@input[[i]])
  load(snakemake@input[[i]])
  txi$counts = calc_rpkm(txi$counts, txi$length)
  save(txi, file = snakemake@output[[i]])
}