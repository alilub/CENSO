#
# Script for plotting barplots showing the library sizes of samples
#
source("scripts/Functions.R")
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library(ggplot2))
#
methods = c("raw", 
            as.vector(read.csv(snakemake@config$methods, sep = "\t")[,1]))
methods.type = c("count", 
                 as.vector(read.csv(snakemake@config$methods, sep = "\t")[,2]))
mode = script_mode[snakemake@params$mode][[1]]
save.image(paste0("images/", snakemake@rule, ".RData"))


out = 1
for(i in 1:length(snakemake@input$raw)){
  #
  # Loading data
  #
  load(snakemake@input$raw[i])
  df.Lib_sizes = colSums(txi$counts)
  df.methode = rep("raw", length(df.Lib_sizes))
  df.methode.kind = rep(methods.type[1], length(df.Lib_sizes))
  for(j in 1:(length(methods)-1)){
    poses = (i-1)*(length(methods)-1)
    load(snakemake@input$normed[j + poses])
    tmp = colSums(txi$counts)
    df.Lib_sizes = c(df.Lib_sizes, tmp)
    tmp = rep(methods[j+1], length(colSums(txi$counts)))
    df.methode = c(df.methode, tmp)
    tmp = rep(methods.type[j+1], length(colSums(txi$counts)))
    df.methode.kind = c(df.methode.kind, tmp)
  }
  #
  # To data.frame
  #
  df.sample = names(df.Lib_sizes)
  species = extend_names(names = df.sample,
                         snakemake = snakemake,
                         parts = mode[[1]])[1]
  df.sample = extend_names(names = df.sample,
                         snakemake = snakemake,
                         parts = mode[[2]])
  df.methode = factor(df.methode, levels = unique(df.methode))
  df.methode.kind = factor(df.methode.kind, levels = unique(df.methode.kind))
  df = data.frame(LibrarySizes = df.Lib_sizes,
                  methode = df.methode,
                  sample = df.sample,
                  kind = df.methode.kind)
  #
  # Summarize plots for each kind of methode
  #
  for(j in sort(unique(as.vector(df.methode.kind)))){
    df.tmp = df[which(df$kind == j),]
    poses = which(methods.type == j)
    bp <- ggplot(data=df.tmp, aes(x=methode, y=LibrarySizes, fill=sample)) +
      geom_bar(stat="identity", position=position_dodge(), color = "black")+ 
      scale_fill_manual(values = get_colour_met(length(unique(df.sample))), name = "sample")+
      theme_linedraw()+
      theme(axis.text.x=element_text(angle = -90, hjust = 0),
            text = element_text(size=20))+
      ggtitle(paste("Library Sizes of", species),
              subtitle = paste("Genes:", length(txi$counts[,1])))
    if(j == "count")
      bp <- bp + geom_hline(yintercept = (10^6), color='red', linetype = "dashed")
    if(snakemake@params$mode == "tissues"){
      png(snakemake@output[[out]], 
          width = 2440,
          res = 120)
      print(bp)
      dev.off()
      out = out + 1
    } else {
      png(snakemake@output[[out]], 
          width = 1440, 
          res = 120)
      print(bp)
      dev.off()
      out = out + 1
    }
  }
}