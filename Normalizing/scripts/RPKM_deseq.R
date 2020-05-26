#
# Script for calculating DESeq normalized RPKM values
#
source("scripts/Functions.R")
save.image(paste0("images/", snakemake@rule, ".RData"))
suppressPackageStartupMessages(library(DESeq2))

for(i in 1:length(snakemake@input)){
  #
  # Plot boxplot of log counts
  #
  print(snakemake@input[[i]])
  load(snakemake@input[[i]])
  RPKM = calc_rpkm(txi$counts, txi$length)
  RPKM = round(RPKM)
  #
  # Convert data structure for DESeq
  #
  sampleTable <- data.frame(condition = extract_cons(colnames(txi$counts)))
  rownames(sampleTable) <- colnames(txi$counts)
  dds <- DESeqDataSetFromMatrix(RPKM, sampleTable, ~condition)
  sf <- estimateSizeFactorsForMatrix(counts(dds))
  #
  # Get normalized counts
  #
  sizeFactors(dds) <- sf
  txi$counts = counts(dds, normalized=TRUE)
  #
  # Saving normalized counts and scale factors
  #
  save(txi, file = snakemake@output[[i]])
  sf = 1/sf
  save(sf, file = snakemake@output[[i+length(snakemake@input)]])
}