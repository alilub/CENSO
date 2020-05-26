#
# Script for creating CoV-Plots summarized for each feat
#

#
# Calculating the CoV values and returning a data frame
#
calc.CoV.all <- function(counts, methode, species, type) {
  cov.methode = vector()
  cov.cov = vector()
  cov.feat = vector()
  cov.species = vector()
  cov.met.type= vector()
  if(snakemake@params$mode == "tissues"){
    samples = get_species(colnames(counts))
  } else {
    samples = get_tissues(colnames(counts))
  }
  for(y in 1:length(counts[,1])){
    cov.feat[length(cov.feat) + 1] = "summary"
    cov.methode[length(cov.methode) + 1] = methode
    cov.cov[length(cov.cov) + 1] = sd(counts[y,])/mean(counts[y,])
    cov.met.type[length(cov.met.type) + 1] = type
    cov.species[length(cov.species) + 1] = species
  }
  cov.cov[is.nan(cov.cov)] = 0
  erg = data.frame(feat = cov.feat,
                   methode = cov.methode,
                   cov = cov.cov,
                   spec = cov.species,
                   type = cov.met.type)
  mean_cov <<- rbind(mean_cov,
                     data.frame(mean_cov = mean(cov.cov), 
                                method = methode,
                                species = species,
                                type = type))
  return(erg)
}

plot_cv <- function(cov, feature, cols){
  cvs = cov.tmp[which(cov$spec == feature),]
  cvs$feat = rep("", length(cvs$feat))
  gg <- ggplot(cvs, aes(feat, cov, fill=methode)) +
    theme_light()+
    geom_boxplot(outlier.shape = NA) + 
    scale_fill_manual(values = cols, name = "methods")+
    ggtitle(feature)+
    xlab("")+
    ylim(c(0,4.5))
  return(gg)
}

wrapper <- function(x, cov, cols){
  if(mode[[1]][1] == 2)
    features = c("Brain", "Cerebellum", "Heart", "Kidney", "Liver", "Testis")
  else
    features = c("Bonobo", "Chicken", "Chimp", "Gorilla", "Human", "Macaque", "Mouse", "Opossum", "Orangutan", "Platypus")
  return(plot_cv(cov = cov, feature = features[x], cols = cols))
}

source("scripts/Functions.R")
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
#
# for species
#
methods = c("raw", 
            as.vector(read.csv(snakemake@config$methods, sep = "\t")[,1]))
methods.type = c("count", 
                 as.vector(read.csv(snakemake@config$methods, sep = "\t")[,2]))
mode = script_mode[snakemake@params$mode][[1]]
save.image(paste0("images/", snakemake@rule, ".RData"))
#
# Calculating
#
cov = data.frame(feat = vector(),
                 methode = vector(),
                 cv = vector(),
                 spec = vector(),
                 type = vector())
mean_cov <<- data.frame(mean_cov = vector(), 
                        method = vector(),
                        species = vector(),
                        type = vector())
#
# Loading and calculating
#
for(i in 1:length(snakemake@input$raw)){
  load(snakemake@input$raw[i])
  species = extend_names(names = colnames(txi$counts),
                         snakemake = snakemake,
                         parts = mode[[1]])[1]
  
  cov = rbind(cov, calc.CoV.all(txi$counts, "raw", species, methods.type[1]))
  for(j in 1:(length(methods)-1)){
    poses = (i-1)*(length(methods)-1)
    load(snakemake@input$normed[j + poses])
    cov = rbind(cov, calc.CoV.all(txi$counts, methods[j+1], species, methods.type[j+1]))
  }
}
cov$methode = factor(cov$methode, levels = unique(cov$methode))
#
# Plotting
#
out = 1
for(j in sort(unique(as.vector(methods.type)))){
  cov.tmp = cov[which(cov$type == j),]
  if(j != "count"){
    cov.tmp.2 = cov[which(cov$methode == "raw"),]
    cov.tmp = rbind(cov.tmp, cov.tmp.2)
    cols = get_colour_met(length(methods))[c(1, which(methods.type == j))]
  } else {
    cols = get_colour_met(length(methods))[which(methods.type == j)]
  }
  
  x = c(1:length(unique(cov.tmp$spec)))
  bp <- ggarrange(plotlist =  lapply(x, wrapper, cols = cols, cov = cov.tmp),
                  common.legend = TRUE, 
                  legend="bottom")
  bp <- annotate_figure(bp,
                        top = text_grob(paste("Coeffienct of Variation")))
  png(snakemake@output[[out]], 
      width = 700,
      height = 1000,
      res = 100)
  print(bp)
  dev.off()
  out = out + 1
}


save(mean_cov, file = snakemake@output[[out]])