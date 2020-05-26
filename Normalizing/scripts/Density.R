#
# Script for creating pca plots for species
#
source("scripts/Functions.R")
#
# libraries
#
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(ggfortify))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library(ggpubr))
#
#
#
methods = c("raw", 
            as.vector(read.csv(snakemake@config$methods, sep = "\t")[,1]))
methods.type = c("count", 
                 as.vector(read.csv(snakemake@config$methods, sep = "\t")[,2]))
mode = script_mode[snakemake@params$mode][[1]]
out = 1
save.image(paste0("images/", snakemake@rule, ".RData"))

dens = list()
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
  
  counts = log(counts[which(!is.na(rownames(counts))),] + 1)
  
  species = extend_names(names = colnames(counts),
                         snakemake = snakemake,
                         parts = c(1))
  tis = extend_names(names = colnames(counts),
                     snakemake = snakemake,
                     parts = c(2))
  sample = extend_names(names = colnames(counts),
                        snakemake = snakemake,
                        parts = c(3))
  
  data = data.frame(value = vector(),
                    species = vector(),
                    tis = vector(),
                    sample = vector())
  for(i in 1:length(species)){
    tmp = data.frame(value = counts[,i],
                     species = species[i],
                     tis = tis[i],
                     sample = sample[i])
    data = rbind(data, tmp)
  }
  
  #
  # Colours
  #
  cols = get_colour_feat(length(unique(sample)))
  
  dens[[length(dens) + 1 ]] = ggplot(data = data, aes(x=value, color = sample)) + 
    geom_density()+
    facet_grid(tis ~species)+
    ggtitle(methods[j])+
    scale_x_continuous(limits = c(-5, 10))+
    theme_linedraw()
}

for (j in sort(unique(sort(methods.type)))){
  poses = which(methods.type == j)
  pcaPlot <- ggarrange(plotlist = dens[poses], 
                       common.legend = TRUE, 
                       legend="bottom")
  pcaPlot <- annotate_figure(pcaPlot,
                             top = text_grob("PCAs of species and tissues"))
  png(snakemake@output[[out]], 
      width = 2880, 
      height = 2880, 
      res = 140)
  print(pcaPlot)
  dev.off()
  out = out + 1
}