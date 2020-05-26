get_species <- function(cons){
  erg = vector()
  for(i in 1:length(cons)){
    tmp = strsplit(cons[i], "_")[[1]]
    erg[i] = tmp[1]
  }
  return(erg)
}

methods = c(as.vector(read.csv(snakemake@config$methods_all, sep = "\t")[,1]))

input = list(all = c(dat = "all_dat",
                     input_eve = "all_in",
                     RData = "all_RData",
                     nwk = "all_nwk",
                     nindiv = "all_nindiv",
                     input = "all",
                     input_raw = "all_raw"))

save.image(paste0("images/", snakemake@rule, ".RData"))

for(x in 1:length(input)){
  for(i in 1:length(snakemake@input[input[[x]]["input_raw"]][[1]])){
    for(j in 1:length(methods)){
      poses = (i-1)*length(methods) - i
      if(j == 1){
        load(snakemake@input[input[[x]]["input_raw"]][[1]][[i]])
        file.in = snakemake@input[input[[x]]["input_raw"]][[1]][[i]]
        out = j + poses + i
      } else {
        load(snakemake@input[input[[x]]["input"]][[1]][[j + poses]])
        file.in = snakemake@input[input[[x]]["input"]][[1]][[j + poses]]
        out = j + poses + i 
      }
      sums = colSums(t(txi$counts))
      sums = which(sums == 0)
      if(length(sums) > 0){
        counts = txi$counts[-sums,]
      } else {
        counts = txi$counts
      }
      
      counts = log(counts + 1)
      #Filtering all zeros
      row_sub = apply(counts, 1, function(row) all(row !=0 ))
      counts = counts[row_sub,]
      
      species.txi = get_species(colnames(counts))
      species.order = c("gga", "mdo", "mmu", "mml", "ggo", "ptr", "ppa", "hsa", "ppy", "oan")
      counts.ord = vector()
      indis = vector()
      for (j in 1:length(species.order)) {
        poses = which(species.order[j] == species.txi)
        indis[j] = length(poses)
        counts.ord = cbind(counts.ord, counts[,poses])
      }
      write.table(nrow(counts), file = snakemake@output[input[[x]]["dat"]][[1]][[out]], quote = F, col.names = F, row.names = F)
      write.table(nrow(counts), file = snakemake@output[input[[x]]["input_eve"]][[1]][[out]], quote = F, col.names = F, row.names = F)
      write.table(counts.ord, file = snakemake@output[input[[x]]["dat"]][[1]][[out]], quote = F, col.names = F, append = T)
      write.table(t(indis), file = snakemake@output[input[[x]]["nindiv"]][[1]][[out]], quote = F, row.names = F, col.names = F)
      file.copy(from = snakemake@config$nwk, to = snakemake@output[input[[x]]["nwk"]][[1]][[out]])
      print(file.in)
      file.copy(from = file.in, to = snakemake@output[input[[x]]["RData"]][[1]][[out]])
      print(snakemake@output[input[[x]]["RData"]][[1]][[out]])
    }
  }
}

  
