#
# Script for plotting RLE-plots
#
source("scripts/Functions.R")
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(EDASeq))
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
  frame.sample = vector()
  frame.methode = vector()
  frame.counts = vector()
  frame.feat = vector()
  frame.type = vector()
  load(snakemake@input$raw[i])
  rle.data = no.plotRLE(txi$counts)
  for(j in 1:(length(methods)-1)){
    poses = (i-1)*(length(methods)-1)
    load(snakemake@input$normed[j + poses])
    rle.data = cbind(rle.data, no.plotRLE(txi$counts))
  }
  species = extend_names(names = colnames(txi$counts),
                         snakemake = snakemake,
                         parts = mode[[1]])[1]
  #
  # rle.data to data.frame
  #
  names = colnames(rle.data)
  rle.data = rle.data[, order(names)]
  names = colnames(rle.data)
  names = extend_names(names = names,
                       snakemake = snakemake,
                       parts = mode[[2]])
  tmp.methode = rep(methods, length(rle.data[1,])/length(methods))
  tmp.methode.type = rep(methods.type, length(rle.data[1,])/length(methods.type))
  for(x in 1:length(rle.data[1,])){
    for (y in 1:length(rle.data[,1])){
      frame.sample[length(frame.sample) + 1] = names[x]
      frame.methode[length(frame.methode) + 1] = tmp.methode[x]
      frame.type[length(frame.type) + 1] = tmp.methode.type[x]
      frame.counts[length(frame.counts) + 1] = rle.data[y,x]
      frame.feat[length(frame.feat) + 1] = strsplit(names[x], " ")[[1]][1]
    }
  }
  df.counts = data.frame(sample = frame.sample, 
                         method = frame.methode, 
                         rle = frame.counts,
                         feat = frame.feat,
                         type = frame.type)
  #
  #
  #
  fa.samples = factor(df.counts$sample)
  df.counts$method = factor(df.counts$method, levels = unique(df.counts$method))
  #
  # merge with cov
  #
  load(snakemake@input$cov)
  mean_cov = mean_cov[which(mean_cov$species == species),]
  tmp = round(mean_cov$mean_cov[match(df.counts$method, mean_cov$method)], 2)
  tmp = paste(df.counts$method, tmp, sep = "-")
  df.counts$method = factor(tmp, levels = unique(tmp))
  #
  # Creating one figure for each kind of methode
  #
  for(x in sort(unique(as.vector(methods.type)))){
    tmp.df.counts <- df.counts[which(df.counts$type == x),]
    suppressWarnings(
    bp <- ggplot(data = tmp.df.counts, aes(method, rle, fill=sample)) +
      geom_boxplot(outlier.shape = NA) + 
      theme_linedraw()+
      scale_fill_manual(values = get_colour_met(length(unique(df.counts$sample))), name = "samples")+
      scale_y_continuous(limits = quantile(df.counts$rle, c(0.1, 0.9)))+ 
      theme(axis.text.x=element_text(angle = -90, hjust = 0))+
      ggtitle(paste("RLE plot for", species), 
              subtitle = paste("Genes:", length(txi$counts[,1])))
    )
    if(snakemake@params$mode == "tissues"){
      png(snakemake@output[[out]], 
          width = 2440,
          res = 120)
      print(bp)
      dev.off()
    } else {
      png(snakemake@output[[out]], 
          width = 1440, 
          res = 120)
      print(bp)
      dev.off()
    }
    out = out + 1
  }
}