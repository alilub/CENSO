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

save.image(paste0("images/", snakemake@rule, ".RData"))

dist_all = dist_all.scal = vector()
for(i in 1:length(snakemake@input$raw)){
    dist = vector()
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
      scount = as.vector(txi$counts)
      dist = cbind(dist, scount)
    }
    dist.scal = (dist*10^6*ncol(txi$counts))/colSums(dist)
    dist_all.scal = rbind(dist_all.scal, dist.scal)
    dist_all = rbind(dist_all, dist)
}

    
dist_all = t(log(dist_all + 1))
rownames(dist_all) = methods
tdist = factoextra::get_dist(dist_all)
png(snakemake@output[[out]])
pheatmap::pheatmap(as.matrix(tdist),
                   main = "Distances of normalization methods")
graphics.off()
out = out + 1

dist_all.scal = t(log(dist_all.scal + 1))
rownames(dist_all.scal) = methods
tdist = factoextra::get_dist(dist_all.scal)
png(snakemake@output[[out]])
pheatmap::pheatmap(as.matrix(tdist),
                   main = "Distances of normalization methods (scaled)")
graphics.off()
