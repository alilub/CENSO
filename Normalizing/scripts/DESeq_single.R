#
# Script for normalizing the given rna-seq-data using the DESeq method
#
source("scripts/Functions.R")
suppressPackageStartupMessages(library(DESeq2))

for(i in 1:length(snakemake@input)){
  load(snakemake@input[[i]])
  #
  # Convert data structure
  #
  sampleTable <- data.frame(condition = extract_cons(colnames(txi$counts)))
  rownames(sampleTable) <- colnames(txi$counts)
  dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
  #
  # Normalization
  #
  nm <- assays(dds)[["avgTxLength"]]
  sf <- estimateSizeFactorsForMatrix(counts(dds) / nm)
  sizeFactors(dds) <- sf
  txi$counts = counts(dds, normalized=TRUE)
  #
  # Saving normalized data and scale factors
  #
  save(txi, file = snakemake@output[[i]])
  sf = 1/sf
  save(sf, file = snakemake@output[[i+length(snakemake@input)]])
}
