#
# Scaling TPM values for multiple species
#


#
# Getting species with maximal annotations
#
con_spec = read.csv(snakemake@config$specs, sep = "\t")
annos = as.vector(con_spec$annoted)
max_annos = which(annos == max(annos))
#
# Getting alpha
#
save.image(paste0("images/", snakemake@rule, ".RData"))
for(i in 1:length(snakemake@input)){
  print(snakemake@input[[i]])
  load(snakemake@input[[i]])
  #
  # Convert data structure
  #
  txi$counts = txi$abundance
  old = txi$counts
  #
  # extract ref
  #
  cn = strsplit(colnames(txi$counts), "_")
  tmp = vector()
  for(j in 1:length(cn)){
    tmp[j] = cn[[j]][1]
  }
  cn = tmp
  ref = which(cn == con_spec$short[max_annos])
  mod = which(cn != con_spec$short[max_annos])
  
  
  alpha = 0
  if(length(ref) == 0){
    alpha = 1
  } else {
    for(j in 1:length(txi$counts[,1])){
    tmp = prod(txi$counts[j,ref])
    alpha = alpha + tmp^(1/length(ref))
    }
    alpha = alpha * 10^(-6)
    write(alpha, file = "alpha_ortho.txt")
  }
  #
  # Modifi tpms
  #
  for(j in 1:length(txi$counts[1,])){
    txi$counts[,j] = txi$counts[,j] * alpha
  }
  save(txi, alpha, file=snakemake@output[[i]])
}
