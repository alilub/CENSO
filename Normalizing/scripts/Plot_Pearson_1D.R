#
# Plotting heatmaps showing the simularity of sample densities
#
source("scripts/Functions.R")
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library(geneplotter))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))
#
#
#
methods = c("raw", 
            as.vector(read.csv(snakemake@config$methods, sep = "\t")[,1]))
methods.type = c("count", 
                 as.vector(read.csv(snakemake@config$methods, sep = "\t")[,2]))
mode = script_mode[snakemake@params$mode][[1]]
save.image(paste0("images/", snakemake@rule, "_", snakemake@wildcards[[1]],  ".RData"))
#
#
#
out = 1
for(i in 1:length(snakemake@input$raw)){
  load(snakemake@input$raw[i])
  # Little bit colour
  cols = get_colour_feat(length(colnames(txi$counts)))
  #
  # Title
  #
  feat = extend_names(names = colnames(txi$counts),
                         snakemake = snakemake,
                         parts = mode[[1]])[1]
  #
  # Plots
  #
  log.counts = to.df(txi$counts, "raw")
  for(j in 1:(length(methods)-1)){
    poses = (i-1)*(length(methods)-1)
    load(snakemake@input$normed[j + poses])
    log.counts = rbind(log.counts, to.df(txi$counts, methods[j+1]))
  }
  
  samples = unique(log.counts$sample)
  
  for(j in 1:length(methods)){
    tmp.counts = log.counts[which(log.counts$methode == methods[j]),]
    
    data = data.frame(pearson = vector(),
                      sample1 = vector(),
                      sample2 = vector())
    
    for(x in 1:length(samples)){
      tmp.counts.1 = tmp.counts[which(tmp.counts$sample == samples[x]),]
      species.1 = strsplit(as.vector(samples)[x], " ")[[1]][1]
      for(y in 1:length(samples)){
        tmp.counts.2 = tmp.counts[which(tmp.counts$sample == samples[y]),]
        species.2 = strsplit(as.vector(samples)[y], " ")[[1]][1]
        if(species.1 == species.2){
          tmp.test = cor.test(tmp.counts.1$logcounts, tmp.counts.2$logcounts, method = snakemake@wildcards[[1]])
          tmp = data.frame(pearson = tmp.test$estimate,
                         sample1 = samples[x],
                         sample2 = samples[y])
          data = rbind(data, tmp)
        }
      }
    }
    
    tmp = unlist(strsplit(as.vector(data$sample1), " "))
    tmp = tmp[seq(from = 1, to = length(tmp), by = 3)]
    data$species = as.factor(tmp)
    
    ggp <- ggplot(data = data, aes(x = species, y = pearson, , fill=species)) +
      geom_boxplot(outlier.colour="black")+
      ggtitle(methods[j])+
      theme_linedraw()+
      theme(axis.text.x=element_text(angle = -90, hjust = 0))+
      scale_fill_manual(values = get_colour_met(length(unique(data$species))), name = "species")+
      ylim(0,1)
    
    if(j == 1){
      file = paste0("data/plots/", snakemake@params$mode, "/", feat, "_raw_pearson_1D.png")
      png(filename = file, 
          width = 720, 
          height = 720, 
          res = 120)
    } else {
      png(snakemake@output[[out]],
          width = 720, 
          height = 720, 
          res = 120)
      out = out + 1
    }
    print(ggp)
    dev.off()
  }
}
