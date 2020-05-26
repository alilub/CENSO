#
# Script for creating pca plots for species
#
source("scripts/Functions.R")
#
# libraries
#
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(ggfortify))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(bioDist))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(ggrepel))
#
#
#

methods = c("raw", 
            as.vector(read.csv(snakemake@config$methods, sep = "\t")[,1]))
methods.type = c("count", 
                 as.vector(read.csv(snakemake@config$methods, sep = "\t")[,2]))
out = 1
save.image(paste0("images/", snakemake@rule, ".RData"))

for(j in 1:length(methods)){
  counts = vector()
  for(i in 1:length(snakemake@input$raw)){
    if(j == 1){
      load(snakemake@input[[i]])
    } else {
      poses = length(snakemake@input$raw) + (i-1)*length(methods) - i
      load(snakemake@input[[j + poses]])
    }
    if(i == 1 || length(row.names(txi$counts)) == length(row.names(counts))){
      counts = cbind(counts, txi$counts)
      next
    }
    if(length(row.names(txi$counts)) > length(row.names(counts))){
      tmp = match(row.names(counts), row.names(txi$counts))
      counts = cbind(counts, txi$counts[tmp,])[which(!is.na(tmp)),]
      next
    }
    if(length(row.names(txi$counts)) < length(row.names(counts))){
      tmp = match(row.names(txi$counts), row.names(counts))
      counts = cbind(counts[tmp,], txi$counts)
      next
    }
  }
  
  counts = counts[which(!is.na(rownames(counts))),]
  scount = t(log(counts + 1))
  species = extend_names(names = rownames(scount),
                         snakemake = snakemake,
                         parts = c(1))
  tis = extend_names(names = rownames(scount),
                     snakemake = snakemake,
                     parts = c(2))
  samples = extend_names(names = rownames(scount),
                         snakemake = snakemake,
                         parts = c(3))
  rownames(scount) = paste(species, tis, samples)
  sample = rownames(scount)
  dataInfo = data.frame(sample, species, tis)
  
  f <- function(xx) nj(spearman.dist(xx))
  tdist = spearman.dist(scount)
  tdist = as.matrix(tdist)
  col.cell = cbind(organism = get_colour_feat(length(unique(species)))[factor(species)],
                   tissues = get_colour_feat(length(unique(tis)))[factor(tis)])
  png(snakemake@output[[out]], 
      width = 720,
      height = 720)
  p <- heatmap.plus::heatmap.plus(tdist,
                             ColSideColors=col.cell,
                             scale="row",
                             margins = c(8,8),
                             col = terrain.colors(256),
                             main = "Heatmap of all samples")
  dev.off()
  out = out + 1
}