#
# Calculating TPM values
#

source("scripts/Functions.R")
save.image(paste0("images/", snakemake@rule, ".RData"))

for(i in 1:length(snakemake@input)){
  load(snakemake@input[[i]])
  txi$counts = calc_tpm(txi$counts, txi$length)
  save(txi, file = snakemake@output[[i]])
}