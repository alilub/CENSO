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
save.image(paste0("images/", snakemake@rule, ".RData"))
out = 1
for(i in 1:length(snakemake@input$raw)){
  pcas = list()
  for(j in 1:length(methods)){
    if(j == 1){
      load(snakemake@input[[i]])
    } else {
      poses = length(snakemake@input$raw) + (i-1)*length(methods) - i
      load(snakemake@input[[j + poses]])
    }
    scount = t(log(txi$counts + 1))
    species = extend_names(names = rownames(scount),
                           snakemake = snakemake,
                           parts = c(1))
    tis = extend_names(names = rownames(scount),
                       snakemake = snakemake,
                       parts = c(2))
    rownames(scount) = extend_names(names = rownames(scount),
                                    snakemake = snakemake,
                                    parts = mode[[2]])
    sample = rownames(scount)
    dataInfo = data.frame(sample, species, tis)
    #
    # Colours
    #
    if(snakemake@params$mode == "species"){
      cols = get_colour_feat(length(unique(tis)))
    } else {
      cols = get_colour_feat(length(unique(species)))
    }
    pcDat <- prcomp(scount)
    
    if(snakemake@params$mode == "species"){
      pcPlot <- autoplot(pcDat,
                         data = dataInfo,
                         colour = "tis",
                         size = 3, frame = TRUE) +
        ggtitle(paste("PCA after using", methods[j]),
                subtitle = paste("Genes:", length(txi$counts[,1])))+
        scale_color_manual(values = cols, name = "tis")
    } else {
      pcPlot <-autoplot(pcDat,
                        data = dataInfo,
                        colour = "species",
                        size = 3, frame = TRUE) +
        ggtitle(paste("PCA after using", methods[j]),
                subtitle = paste("Genes:", length(txi$counts[,1])))+
        scale_color_manual(values = cols, name = "species")
    }
    pcas[[length(pcas) + 1 ]] = pcPlot +
      theme_linedraw()
  }
  
  for (j in sort(unique(as.vector(methods.type)))){
    poses = which(methods.type == j)
    pcaPlot <- ggarrange(plotlist = pcas[poses], 
                         common.legend = TRUE, 
                         legend="bottom")
    if(snakemake@params$mode == "species"){
      pcaPlot <- annotate_figure(pcaPlot,
                                 top = text_grob(paste("PCAs of", species[1])))
    } else {
      pcaPlot <- annotate_figure(pcaPlot,
                                 top = text_grob(paste("PCAs of", tis[1])))
    }
    png(snakemake@output[[out]], 
        width = 1440, 
        height = 1440, 
        res = 140)
    print(pcaPlot)
    dev.off()
    out = out + 1
  }
}