#
# Script for creating a data set containing all orthologs genes from given species as one txi object
#
orthos = read.csv(snakemake@input$orthos, sep = ",")
#
# Build list
#
txi_new = list(matrix(nrow = length(orthos[,1]), ncol = 0),
               matrix(nrow = length(orthos[,1]), ncol = 0),
               list(),
               matrix(nrow = length(orthos[,1]), ncol = 0),
               "no")
names(txi_new) = c("abundance", "counts", "infReps", "length", "countsFromAbundance")
miss = vector()
for(i in 1:length(snakemake@input$reads)){
  print(snakemake@input[[i]])
  load(snakemake@input[[i]])
  #
  # get right collum
  #
  genes = row.names(txi$counts)
  j = 0
  poses = vector()
  for(j in 1:length(orthos[1,])){
    poses = match(orthos[,j], genes)
    tmp = poses[!is.na(poses)]
    if(length(tmp) > 0)
      break
  }
  miss = c(miss, which(is.na(poses)))
  #
  # matching
  #
  # abundance
  tmp = txi$abundance[poses,]
  row.names(tmp) = orthos[,1]
  tmp[is.na(tmp)] = 0
  txi_new$abundance = cbind(txi_new$abundance, tmp)
  # counts
  tmp = txi$counts[poses,]
  row.names(tmp) = orthos[,1]
  tmp[is.na(tmp)] = 0
  txi_new$counts = cbind(txi_new$counts, tmp)
  # length
  tmp = txi$length[poses,]
  row.names(tmp) = orthos[,1]
  tmp[is.na(tmp)] = 0
  txi_new$length = cbind(txi_new$length, tmp)
}
#
# Deleting missing data
#
miss = unique(miss)
if(length(miss) > 0){
  # abundance
  txi_new$abundance = txi_new$abundance[-miss,]
  # counts
  txi_new$counts = txi_new$counts[-miss,]
  # length
  txi_new$length = txi_new$length[-miss,]
}
txi = txi_new
remove(txi_new)
# Saving
save(txi, file = snakemake@output[[1]])