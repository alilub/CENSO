#
# Splitting a txi object containing all count data into tissues
#
source("scripts/Functions.R")
tissues = read.csv(file = snakemake@config$tissues, sep = "\t")
load(snakemake@input[[1]])
txi_all = txi
mode = snakemake@params$mode
save.image(paste0("images/", snakemake@rule, ".RData"))

for(i in 1:length(tissues[,1])){
  print(as.vector(tissues[i,1]))
  tissue.short = as.vector(tissues[i,2])
  samples = colnames(txi_all$counts)
  samples.tissues = get_tissues(samples)
  poses = which(samples.tissues == tissue.short)
  txi = list(abundance = txi_all$abundance[,poses],
             counts = txi_all$counts[,poses],
             infReps = txi_all$infReps,
             length = txi_all$length[,poses],
             countsFromAbundance = txi_all$countsFromAbundance)
  #
  # Filtering missing data
  #
  if(mode != "all"){
    count = t(txi$counts)
    species = extend_names(names = colnames(txi$counts), 
                           snakemake = snakemake, 
                           parts = c(1))
    species = factor(species, levels = unique(species))
    sums = rowsum(count, group = species)
    poses = vector()
    for(x in 1:length(sums[1,])){
      if(sum(sums[,x]) == 0)
        next
      tmp = which(sums[,x] == 0)
      if((length(tmp)/2) <= length(sums[,x]))
        poses[length(poses) + 1] = x
    }
    #
    # saving
    #
    txi = list(abundance = txi$abundance[poses,],
               counts = txi$counts[poses,],
               infReps = txi$infReps,
               length = txi$length[poses,],
               countsFromAbundance = txi$countsFromAbundance)
  }
  
  save(txi, file = snakemake@output[[i]])
}