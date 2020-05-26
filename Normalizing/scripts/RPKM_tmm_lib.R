#
# Calculating RPKM values with TMM-normalized library sizes
#
source("scripts/Functions.R")
suppressPackageStartupMessages(library(edgeR))

save.image(paste0("images/", snakemake@rule, ".RData"))

for(i in 1:length(snakemake@input)){
  load(snakemake@input[[i]])
  cts <- txi$counts
  normMat <- txi$length
  normMat <- normMat/exp(rowMeans(log(normMat)))
  o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
  groups = extract_cons(colnames(txi$counts))
  #
  # One per species
  #
  species = get_species(colnames(txi$counts))
  s.species = unique(species)
  RPKM = matrix(ncol = 0, nrow = length(txi$counts[,1]))
  rownames(RPKM) = rownames(txi$counts)
  for(x in 1:length(s.species)){
    poses = which(species == s.species[x])
    y <- DGEList(cts[,poses], group = factor(groups[poses]))
    y$offset <- t(t(log(normMat)) + o)[,poses]
    y = calcNormFactors(y, method = "TMM")
    RPKM = cbind(RPKM, rpkm(y, gene.length = txi$length[,poses[1]]))
  }
  names = colnames(txi$counts)
  txi$counts = RPKM
  colnames(txi$counts) = names
  save(txi, file = snakemake@output[[i]])
}