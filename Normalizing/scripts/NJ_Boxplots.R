#
# Function for plotting phylogenetic tree using spearman correlation as distance
# For species
#
sep_feat = function(s){
  tmp = strsplit(s, split = "/", fixed = T)[[1]][4]
  tmp = strsplit(tmp, split = "_", fixed = t)[[1]][1]
  return(tmp)
}

sep_feat_name = function(s){
  erg = vector()
  tmp = strsplit(s, split = " ")
  for (i in 1:length(tmp)) {
    erg = c(erg, tmp[[i]][1])
  }
  return(erg)
}

bp_tl <- function(data1, data2){
  erg = ggplot(data = data1, aes(feat, length, col = feat)) +
    geom_boxplot(outlier.shape = NA)+
    theme(axis.text.x=element_text(angle = -90, hjust = 0))+
    scale_fill_manual(values = get_colour_feat(length(data1$feat)), name = "feat")+
    ggtitle(paste("Tree Length of", data1$method[1]))+
    theme_linedraw() +
    geom_point(data = data2, aes(feat.o, length, col = feat.o))
  return(erg)
}

wrapper <- function(x){
  return(bp_tl(tree.length[which(tree.length$method == methods[x]), ],
               tree.length.o[which(tree.length.o$method == methods[x]), ])
         )
}

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
tree.length.o = data.frame(length = vector(),
                           feat.o = vector(),
                           method = vector())
tree.length = data.frame(length = vector(),
                         feat = vector(),
                         method = vector())
for(i in 1:length(snakemake@input)){
  #
  # Loading data
  #
  load(snakemake@input[[i]])
  for(x in 1:length(list.tree)){
    depth = ape::node.depth.edgelength(list.tree[[x]])
    tmp = data.frame(length = max(depth),
                     feat.o = sep_feat(snakemake@input[[i]]),
                     method = names(list.tree)[x])
    tree.length.o = rbind(tree.length.o, tmp)
  }
  
  for(x in 1:length(list.bp.tree)){
    for(y in 1:length(list.bp.tree[[x]]$trees)){
      depth = ape::node.depth.edgelength(list.bp.tree[[x]]$trees[[y]])
      tmp = data.frame(length = max(depth),
                     feat = sep_feat(snakemake@input[[i]]),
                     method = names(list.bp.tree)[x])
      tree.length = rbind(tree.length, tmp)
    }
  }
}

for(i in 1:length(methods)){
  plots = ggarrange(plotlist = lapply(i, wrapper), common.legend = T)
  png(snakemake@output[[out]], height = 700, width = 700)
  print(plots)
  dev.off()
  out = out + 1
}
