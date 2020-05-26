#
# Script for deleting genes without reads in any sample
#
source("scripts/Functions.R")

for(i in 1:length(snakemake@input)){
  load(snakemake@input[[i]])
  #
  # delete those genes without any counts
  #
  tcounts = t(txi$counts)
  sums = colSums(tcounts)
  misses = which(sums == 0)
  if(length(misses) > 0){
    # abundance
    txi$abundance = txi$abundance[-misses,]
    # counts
    txi$counts = txi$counts[-misses,]
    # length
    txi$length = txi$length[-misses,]
    # tmp
    tcounts = tcounts[,-misses]
  }
  #
  # Correction of TPM values
  #
  txi$abundance = calc_tpm(txi$counts, txi$length)
  #
  # Saving
  #
  save(txi, file = snakemake@output[[i]])
}