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
suppressPackageStartupMessages(library(ggrepel))
#
#
#
methods = c("raw", 
            as.vector(read.csv(snakemake@config$methods, sep = "\t")[,1]))
methods.type = c("count", 
                 as.vector(read.csv(snakemake@config$methods, sep = "\t")[,2]))
mode = script_mode[snakemake@params$mode][[1]]
out = 1
out.txt = 1

save.image(paste0("images/", snakemake@rule, ".RData"))

list.tree <- list()
list.labels <- list()
list.cols <- list()
list.bp.tree <- list()

plot_dendogram <- function(tree, labels, cols, method){
  p <- ggtree(tree, branch.length = "none") +
    geom_label_repel(aes(label = labels, fill = set_label_levels(label)), colour = cols)+
    scale_fill_manual(values = c("grey90", "white", "yellow", "orange"))+
    ggtitle(method)+
    theme_tree2() +
    theme(legend.position = "none") +
    scale_y_reverse()
  return(p)
}

plot_tree <- function(tree, labels, cols, method){
  p <- ggtree(tree) +
    geom_label_repel(aes(label = labels, fill = set_label_levels(label)), colour = cols) +
    scale_fill_manual(values = c("grey90", "white", "yellow", "orange")) +
    ggtitle(method) +
    theme_tree2() +
    theme(legend.position = "none") +
    scale_y_reverse()
  return(p)
}

wrapper <- function(x){
  if(snakemake@params$plot == "dendo")
    return(plot_dendogram(list.tree[[x]], list.labels[[x]], list.cols[[x]], methods[x]))
  if(snakemake@params$plot == "tree")
    return(plot_tree(list.tree[[x]], list.labels[[x]], list.cols[[x]], methods[x]))
}

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
    
    tree = NULL
    if(snakemake@params$dis == "eukl"){
      f <- function(xx) nj(get_dist(xx))
      tdist = factoextra::get_dist(dat)
    }
    if(snakemake@params$dis == "spearman"){
      f <- function(xx) nj(spearman.dist(xx))
      tdist = spearman.dist(dat)
    }
    if(snakemake@params$cluster == "hclust"){
      tree = hclust(tdist)
      tree = as.phylo(tree)
    }
    if(snakemake@params$cluster == "nj"){
      tree = nj(tdist)
    }
    if(is.null(tree)){
      print("Wrong parameter submitted!")
      q()
    }
    
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
      test = try(root(tree, outgroup))
      if(class(test) == "phylo"){
        tree = test
      }
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
    
    list.tree[[j]] <- tree
    list.labels[[j]] <- labels
    list.cols[[j]] <- cols
    list.bp.tree[[j]] <- bp_trees
  }
  
  feat_2 = extend_names(names = colnames(txi$counts),
                        snakemake = snakemake,
                        parts = mode[[1]])
  
  names(list.tree) = methods
  names(list.bp.tree) = methods
  save(list.tree, list.bp.tree, file = snakemake@output$data[[i]])
  
  for(j in 1:length(sort(unique(as.vector(methods))))){
    #poses = which(methods.type == unique(methods.type)[j])
    pltrees <- ggarrange(plotlist = lapply(j, wrapper))
    pltrees <- annotate_figure(pltrees,
                               top = text_grob(paste("Phylograms of", feat_2)))
    print(snakemake@output$plots[[out]])
    print(out)
    png(snakemake@output$plots[[out]], 
        width = 1440, 
        height = 1440, 
        res = 120)
    print(pltrees)
    dev.off()
    out = out + 1
  }
  
  sink(snakemake@output$text[[out.txt]], append=FALSE, split=FALSE)
  for(j in 1:length(sort(unique(as.vector(methods))))){
    for(k in j:length(sort(unique(as.vector(methods))))){
      print(methods[j])
      print(methods[k])
      print(all.equal(list.tree[[j]], list.tree[[k]]))
      print(all.equal.list(list.tree[[j]], list.tree[[k]]))
    }
  }
  out.txt = out.txt + 1
}