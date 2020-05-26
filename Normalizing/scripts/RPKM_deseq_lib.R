#
# Script for calculating RPKM values using DESeq normalized Library sizes
#
source("scripts/Functions.R")
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(edgeR))
save.image(paste0("images/", snakemake@rule, ".RData"))


for(i in 1:length(snakemake@input)){
    print(snakemake@input[[i]])
    load(snakemake@input[[i]])
    #
    # Convert data structure for DESeq
    #
    sampleTable <- data.frame(condition = extract_cons(colnames(txi$counts)))
    rownames(sampleTable) <- colnames(txi$counts)
    dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
    nm <- assays(dds)[["avgTxLength"]]
    sf <- estimateSizeFactorsForMatrix(counts(dds) / nm)
    #
    # Convert do DGEList
    #
    cts <- txi$counts
    normMat <- txi$length
    normMat <- normMat/exp(rowMeans(log(normMat)))
    o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
    groups = extract_cons(colnames(txi$counts))
    #
    # RPKM normalization
    #
    species = get_species(colnames(txi$counts))
    s.species = unique(species)
    RPKM = matrix(ncol = 0, nrow = length(txi$counts[,1]))
    rownames(RPKM) = rownames(txi$counts)
    for(x in 1:length(s.species)){
      poses = which(species == s.species[x])
      y <- DGEList(cts[,poses], group = factor(groups[poses]), norm.factors = sf[poses])
      y$offset <- t(t(log(normMat)) + o)[,poses]
      RPKM = cbind(RPKM, rpkm(y, gene.length = txi$length[,poses[1]]))
    }
    names = colnames(txi$counts)
    txi$counts = RPKM
    #
    # saving
    #
    txi$counts = fpkm(dds)
    save(txi, file = snakemake@output[[i]])
}
