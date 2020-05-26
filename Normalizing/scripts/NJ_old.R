#
# Function for plotting phylogenetic tree using spearman correlation as distance
# For species
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
  #
  # Plotting
  #
  for(j in 1:length(methods)){
    if(j == 1){
      load(snakemake@input[[i]])
    } else {
      poses = length(snakemake@input$raw) + (i-1)*length(methods) - i
      print(snakemake@input[[j + poses]])
      load(snakemake@input[[j + poses]])
    }
    scount = txi$counts
    colnames(scount) = extend_names(names = colnames(scount),
                                    snakemake = snakemake,
                                    parts = mode[[2]])
    dat = log(t(scount) + 1)
    f <- function(xx) nj(spearman.dist(xx))
    tdist = spearman.dist(dat)
    tree = nj(tdist)
    
    bp <- boot.phylo(phy =  tree, 
                     x = dat, 
                     FUN =  f, 
                     quiet = T, 
                     B = 1000)
    tree$node.label = (bp/1000)*100
    
    if(snakemake@params$mode == "species"){
      feat = extend_names(names = colnames(txi$counts),
                        snakemake = snakemake,
                        parts = c(2))
      cols = get_colour_feat(length(unique(feat)))[factor(feat)]
    } else {
      outgroup = colnames(scount)[grep("Chicken", colnames(scount))]
      tree= root(tree, outgroup)
      #
      # A little bit of colour
      #
      species = extend_names(names = colnames(txi$counts),
                             snakemake = snakemake,
                             parts = c(1))
      cols = get_colour_feat(length(unique(species)))[factor(species)]
    }
    
    bp_trees <- boot.phylo(phy =  tree,
                           x = dat,
                           FUN =  f,
                           quiet = T,
                           B = 1000,
                           trees = T)
    
    cols = c(cols[match(tree$tip.label, colnames(scount))],
             rep("black", length(tree$node.label)))
    labels = c(colnames(scount), tree$node.label)
    
    plot_tree(tree, labels, cols, bp_trees, methods[j])
  }
  
  feat_2 = extend_names(names = colnames(txi$counts),
                         snakemake = snakemake,
                         parts = mode[[1]])
  
  save(list.tree, list.bp.tree, file = snakemake@output$data[[i]])
  
  for(j in 1:length(unique(methods.type))){
    poses = which(methods.type == unique(methods.type)[j])
    pltrees <- ggarrange(plotlist = list.plots[poses])
    pltrees <- annotate_figure(pltrees,
                               top = text_grob(paste("Phylograms of", feat_2)))
    png(snakemake@output$plots[[out]], 
        width = 2880, 
        height = 2880, 
        res = 100)
    print(pltrees)
    dev.off()
    out = out + 1
  }
}