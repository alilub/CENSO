get_species <- function(cons){
  erg = vector()
  for(i in 1:length(cons)){
    tmp = strsplit(cons[i], "_")[[1]]
    erg[i] = tmp[1]
  }
  return(erg)
}
dirs = list.dirs(recursive = F)

for(d in 1:length(dirs)){
  fil = list.files(path = dirs[d], pattern = "RData", full.names = T)
  load(fil)
  sums = colSums(t(txi$counts))
  sums = which(sums == 0)
  if(length(sums) > 0){
    counts = txi$counts[!sums,]
  } else {
    counts = txi$counts
  }
  
  species.txi = get_species(colnames(counts))
  species.order = c("gga", "mdo", "mmu", "mml", "ggo", "ptr", "ppa", "hsa", "ppy", "oan")
  counts.ord = vector()
  indis = vector()
  for (j in 1:length(species.order)) {
    poses = which(species.order[j] == species.txi)
    indis[j] = length(poses)
    counts.ord = cbind(counts.ord, counts[,poses])
  }
  write.table(nrow(counts), file = paste(dirs[d], "Expressions.dat", sep = "/"), quote = F, col.names = F, row.names = F)
  write.table(counts.ord, file = paste(dirs[d], "Expressions.dat", sep = "/"), quote = F, col.names = F, append = T)
  write.table(t(indis), file = paste(dirs[d], "Species.nindiv", sep = "/"), quote = F, row.names = F, col.names = F)
  file.copy(from = "../Species.nwk", to = paste(dirs[d], "Species.nwk", sep = "/"))
}


