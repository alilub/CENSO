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
save.image(paste0("images/", snakemake@rule, ".RData"))
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
    
    data = data.frame(stat = vector(),
                      p.value = vector(),
                      sample1 = vector(),
                      sample2 = vector())
    for(x in 1:length(samples)){
      tmp.counts.1 = tmp.counts[which(tmp.counts$sample == samples[x]),]
      for(y in 1:length(samples)){
        tmp.counts.2 = tmp.counts[which(tmp.counts$sample == samples[y]),]
        ks = ks.test(tmp.counts.1$logcounts, tmp.counts.2$logcounts)
        wm = wilcox.test(tmp.counts.1$logcounts, tmp.counts.2$logcounts)
        if(ks$p.value < 0.05)
          p = "*"
        else
          p = ""
        tmp = data.frame(stat = ks$statistic,
                         p.value = p,
                         sample1 = samples[x],
                         sample2 = samples[y])
        data = rbind(data, tmp)
      }
    }
    
    ggp <- ggplot(data = data, aes(x = sample1, y = sample2, , fill=stat)) +
      geom_tile()+
      ggtitle(methods[j])+
      geom_text(aes(label = p.value))+
      theme(axis.text.x=element_text(angle = -90, hjust = 0))+
      scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                           midpoint = 0.5, limit = c(0,1), 
                           space = "Lab", 
                           name="Two sides\nKolmogorov-mirnov test")
    
    if(j == 1){
      file = paste0("data/plots/", snakemake@params$mode, "/", feat, "_raw_ks.png")
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
