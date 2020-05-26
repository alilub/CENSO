#
# Script for creating dendograms with euklidian or spearman distance (depend on given parameter)
# For species
#
source("scripts/Functions.R")
#
# libraries
#
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(ggfortify))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggdendro))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(bioDist))
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
# Variable for maniging the output
#
out = 1
for(i in 1:length(snakemake@input$raw)){
  #
  # Creating one dendrogram per methode
  #
  for(j in 1:length(methods)){
    #
    # Loading data
    #
    if(j == 1){
      load(snakemake@input[[i]])
    } else {
      poses = length(snakemake@input$raw) + (i-1)*length(methods) - i
      load(snakemake@input[[j + poses]])
    }
    scount = txi$counts
    colnames(scount) = extend_names(names = colnames(scount),
                                    snakemake = snakemake,
                                    parts = mode[[2]])
    #
    # Calculating distance matrix
    #
    tree = NULL
    if(snakemake@params$dis == "eukl"){
      tree = get_dist(log(t(scount) + 1))
    }
    if(snakemake@params$dis == "spearman"){
      tree = spearman.dist(log(t(scount) + 1))
    }
    if(is.null(tree)){
      print("Wrong parameter submitted!")
      q()
    }
    #
    # Clustering
    #
    tree = hclust(tree)
    tree = dendro_data(tree)
    feat = extend_names(names = colnames(txi$counts),
                        snakemake = snakemake,
                        parts = mode[[2]][1])
    cols = get_colour_feat(length(unique(feat)))[factor(feat)]
    cols = cols[match(tree$labels$label, colnames(scount))]
    
    #
    # Plotting dendogram and saving the plot
    #
    plot_cluster(tree = tree, 
                 cols = cols, 
                 method = methods[j])
  }
  
  feat_2 = extend_names(names = colnames(txi$counts),
                         snakemake = snakemake,
                         parts = mode[[1]])
  #
  # Collect all dendograms for one kind of methode an printing the plots
  #
  for (j in sort(unique(as.vector(methods.type)))){
    poses = which(methods.type == j)
    pltrees <- ggarrange(plotlist = list.plots[poses])
    pltrees <- annotate_figure(pltrees,
                               top = text_grob(paste("Dendrograms of", feat_2)))
    png(snakemake@output[[out]], 
        width = 2400, 
        height = 1440, 
        res = 120)
    print(pltrees)
    dev.off()
    out = out + 1
  }
}