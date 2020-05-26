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

norm.vec <- function(x) sqrt(sum(x^2))
print_plot <- function(plot, file){
  png(file, 
      width = 1440, 
      height = 1440, 
      res = 120)
  print(plot)
  dev.off()
}
cluster.variance = data.frame(cluster = vector(),
                              method = vector(),
                              inter = vector(),
                              intra = vector())

methods = c("raw", 
            as.vector(read.csv(snakemake@config$methods, sep = "\t")[,1]))
methods.type = c("count", 
                 as.vector(read.csv(snakemake@config$methods, sep = "\t")[,2]))
mode = script_mode[snakemake@params$mode][[1]]
out = 1
save.image(paste0("images/", snakemake@rule, ".RData"))

pcas = list()
for(j in 1:length(methods)){
  counts = vector()
  for(i in 1:length(snakemake@input$raw)){
    if(j == 1){
      load(snakemake@input[[i]])
    } else {
      poses = length(snakemake@input$raw) + (i-1)*length(methods) - i
      load(snakemake@input[[j + poses]])
    }
    if(i == 1 || length(row.names(txi$counts)) == length(row.names(counts))){
      counts = cbind(counts, txi$counts)
      next
    }
    if(length(row.names(txi$counts)) > length(row.names(counts))){
      tmp = match(row.names(counts), row.names(txi$counts))
      counts = cbind(counts, txi$counts[tmp,])[which(!is.na(tmp)),]
      next
    }
    if(length(row.names(txi$counts)) < length(row.names(counts))){
      tmp = match(row.names(txi$counts), row.names(counts))
      counts = cbind(counts[tmp,], txi$counts)
      next
    }
  }
  
  counts = counts[which(!is.na(rownames(counts))),]
  scount = t(log(counts + 1))
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
  # Calculation of intra und inter cluster variance
  #
  norms = matrix(data = NA, nrow = ncol(counts), ncol = ncol(counts))
  colnames(norms) = colnames(counts)
  rownames(norms) = colnames(counts)
  
  for(x in 1:ncol(counts)){
    for(y in 1:ncol(counts)){
      norms[x,y] = norm.vec(log(counts[,x] + 1) -  log(counts[,y] + 1))
    }
  }
  
  for(x in 1:length(unique(tis))){
    pos = which(tis == unique(tis)[x])
    tmp = data.frame(cluster = "tissue",
                     method = methods[j],
                     intra = sum(norms[pos,pos]^2),
                     inter = sum(norms[-pos,-pos]^2))
    cluster.variance = rbind(cluster.variance, tmp)
  }
  for(x in 1:length(unique(species))){
    pos = which(species == unique(species)[x])
    tmp = data.frame(cluster = "species",
                     method = methods[j],
                     intra = sum(norms[pos,pos]^2),
                     inter = sum(norms[-pos,-pos]^2))
    cluster.variance = rbind(cluster.variance, tmp)
  }
  #
  # Colours
  #
  if(snakemake@params$mode == "species"){
    cols = get_colour_feat(length(unique(tis)))
  } else {
    cols = get_colour_feat(length(unique(species)))
  }
  pcDat <- prcomp(scount)
  
  pcPlot <- autoplot(pcDat,
                     data = dataInfo,
                     colour = "tis",
                     shape = "species",
                     size = 3, frame = TRUE) +
    ggtitle(paste("PCA after using", methods[j]),
            subtitle = paste("Genes:", length(scount[1,])))+
    scale_color_manual(values = cols, name = "tis")+
    scale_shape_manual(values = c(4, 16, 3, 1, 0, 6, 17, 15, 2, 18))+
    theme_linedraw()
  pcas[[length(pcas) + 1 ]] = pcPlot
}

for (j in sort(unique(as.vector(methods.type)))){
  poses = which(methods.type == j)
  pcaPlot <- ggarrange(plotlist = pcas[poses], 
                       common.legend = TRUE, 
                       legend="bottom")
  pcaPlot <- annotate_figure(pcaPlot,
                             top = text_grob("PCAs of species and tissues"))
  print_plot(plot = pcaPlot, file = snakemake@output[[out]])
  out = out + 1
}

pcaPlot <- ggplot(data=cluster.variance, aes(x=method, y=inter, fill=method)) +
  geom_boxplot()+
  scale_fill_manual(values = get_colour_met(length(methods)), name = "methods")+
  theme_linedraw()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0),
        text = element_text(size=10))+
  ggtitle(paste("Inter Cluster variances"))+
  facet_grid(cluster ~.)

print_plot(plot = pcaPlot, file = snakemake@output[[out]])
out = out + 1

pcaPlot <- ggplot(data=cluster.variance, aes(x=method, y = intra, fill=method)) +
  geom_boxplot()+
  scale_fill_manual(values = get_colour_met(length(methods)), name = "methods")+
  theme_linedraw()+
  theme(axis.text.x=element_text(angle = -90, hjust = 0),
        text = element_text(size=10))+
  ggtitle(paste("Tntra Cluster variances"))+
  facet_grid(cluster ~.)

print_plot(plot = pcaPlot, file = snakemake@output[[out]])
out = out + 1