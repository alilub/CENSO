#
# Script for creating CoV-Plots for each sample
#

#
# Calculating the CoV values and returning a data frame
#
calc.CoV <- function(counts, methode, type) {
  if(snakemake@params$mode == "tissues"){
    config = read.csv(snakemake@config$specs, sep = "\t")
    short.tis = as.vector(config$short)
    long.tis = as.vector(config$spec)
  } else
  {
    config = read.csv(snakemake@config$tissues, sep = "\t")
    short.tis = as.vector(config$short)
    long.tis = as.vector(config$tissues)
  }
  cov.feat = vector()
  cov.methode = vector()
  cov.cov = vector()
  cov.met.type= vector()
  for(x in 1:length(short.tis)){
    if(snakemake@params$mode == "tissues"){
      samples = get_species(colnames(counts))
    } else {
      samples = get_tissues(colnames(counts))
    }
    pos = which(samples == short.tis[x])
    if(length(pos) <= 1)
      next
    counts.tmp = counts[,pos]
    for(y in 1:length(counts.tmp[,1])){
      cov.feat[length(cov.feat) + 1] = long.tis[x]
      cov.methode[length(cov.methode) + 1] = methode
      cov.met.type[length(cov.met.type) + 1] = type
      cov.cov[length(cov.cov) + 1] = sd(counts.tmp[y,])/mean(counts.tmp[y,])
    }
  }
  cov.cov[is.nan(cov.cov)] = 0
  erg = data.frame(feat = cov.feat,
                  methode = cov.methode,
                  type = cov.met.type,
                  cov = cov.cov)
  return(erg)
}

source("scripts/Functions.R")
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library(ggplot2))
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
  #
  # Loading and calculating
  #
  load(snakemake@input$raw[i])
  cov = calc.CoV(txi$counts, "raw", "count")
  
  for(j in 1:(length(methods)-1)){
    poses = (i-1)*(length(methods)-1)
    load(snakemake@input$normed[j + poses])
    cov = rbind(cov, calc.CoV(txi$counts, methods[j+1], methods.type[j+1]))
  }
  
  species = extend_names(names = colnames(txi$counts),
                         snakemake = snakemake,
                         parts = mode[[1]])[1]
  cov$methode = factor(cov$methode, levels = unique(cov$methode))
  cov$type= factor(cov$type, levels = unique(cov$type))
  
  for(j in sort(unique(as.vector(cov$type)))){
    cov.tmp = cov[which(cov$type == j),]
    if(j != "count"){
      cov.tmp.2 = cov[which(cov$methode == "raw"),]
      cov.tmp = rbind(cov.tmp, cov.tmp.2)
      poses = c(1, which(methods.type == j))
    } else {
      poses = which(methods.type == j)
    }
    
    bp <- ggplot(data = cov.tmp, aes(feat, cov, fill=methode)) +
      theme_linedraw()+
      geom_boxplot(outlier.shape = NA) + 
      scale_fill_manual(values = get_colour_met(length(methods))[poses], name = "methods")+
      theme(axis.text.x=element_text(angle = -90, hjust = 0))+
      ggtitle(paste("Coefficient of Variantion", species), 
              subtitle = paste("Genes:", length(txi$counts[,1])))
    #
    # Output
    #
    print(snakemake@output[[out]])
    png(snakemake@output[[out]], 
        width = 720, 
        res = 100)
    print(bp)
    dev.off()
    out = out + 1
  }
}