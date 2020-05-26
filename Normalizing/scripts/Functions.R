#
# Collection of functions and variablesused by several scripts
#

#
# Variables
#
script_mode = list(species = list(c(1), c(2,3)),
                   orthos = list(c(1), c(2,3)),
                   tissues = list(c(2), c(1,3)),
                   all = list(c(2), c(1,3)))

#
# Methods
#

#
# For calculating RPKM values
#
calc_rpkm <- function(counts, gens.len){
  #
  # Prepare
  #
  groups = extract_cons(colnames(counts))
  species = get_species(colnames(counts))
  s.species = unique(species)
  RPKM = matrix(ncol = 0, nrow = length(counts[,1]))
  rownames(RPKM) = rownames(counts)
  #
  # Calculate for each species
  #
  for(x in 1:length(s.species)){
    poses = which(species == s.species[x])
    y <- edgeR::DGEList(counts[,poses], group = factor(groups[poses]))
    RPKM = cbind(RPKM, edgeR::rpkm(y, gene.length = gens.len[,poses[1]], normalized.lib.sizes = F))
  }
  #
  # Finishing
  #
  colnames(RPKM) = colnames(counts)
  return(RPKM)
}

#
# function written by Alexander Gabel
#
# For calculating and reaturning tmm-normalized data
calc_tmm <- function(y, saving = T){
  
  tmm_scaleFactors <- y$samples$lib.size * y$samples$norm.factors
  tmm_normFactors <- tmm_scaleFactors/exp(mean(log(tmm_scaleFactors)))
  
  counts_tmm <- as.matrix(y$counts) %*% diag(1/tmm_normFactors)
  colnames(counts_tmm) <- colnames(y$count)
  
  counts_tmm <- round(counts_tmm)
  #
  # Litte extension
  #
  if(saving){
    sf = 1/tmm_normFactors
    names(sf) = colnames(txi$counts)
    save(sf, file = snakemake@output[[i+length(snakemake@input)]])
  }
  #
  # End of extension
  #
  return(counts_tmm)
}

#
# For calculating TPM values
#
calc_tpm <- function(counts, gens.len){
  names = colnames(counts)
  c = counts / gens.len
  TPM = t(t(c)/colSums(c)) * (10^6)
  colnames(TPM) = names
  return(TPM)
}

#
# Returns the long version for sample names given as input, including the parts specifed
#
extend_names <- function(names, snakemake, parts) {
  names_spec = get_species(names)
  names_tis = get_tissues(names)
  names_rep = get_rep(names)
  species = read.csv(snakemake@config$specs, sep = "\t")
  tissues = read.csv(snakemake@config$tissues, sep = "\t")
  for(k in 1:length(names)){
    names_spec[k] = as.vector(species$spec)[which(species$short == names_spec[k])]
    names_tis[k] = as.vector(tissues$tissues)[which(tissues$short == names_tis[k])]
  }
  names = cbind(names_spec, names_tis, names_rep)
  erg = vector()
  for(k in 1:length(names_spec)){
    erg[k] = paste(as.vector(names[k,parts]), collapse=" ")
  }
  return(erg)
}

#
# Returns part of sample names (species and tissue)
#
extract_cons <- function(cons){
  erg = vector()
  for(i in 1:length(cons)){
    tmp = strsplit(cons[i], "_")[[1]]
    erg[i] = paste(tmp[1], tmp[2], sep = "_")
  }
  return(erg)
}

#
# Return vector with colour values for species or tissues for plotting
#
get_colour_feat <- function(len){
  suppressPackageStartupMessages(library(RColorBrewer))
  colourCount = len
  getPalette = colorRampPalette(brewer.pal(12, "Paired"))
  return(getPalette(colourCount))
}

#
# Return vector with colour values for methods for plotting
#
get_colour_met <- function(len){
  suppressPackageStartupMessages(library(RColorBrewer))
  colourCount = len
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  return(getPalette(colourCount))
}

#
# Returning the replicate
#
get_rep <- function(cons){
  erg = vector()
  for(i in 1:length(cons)){
    tmp = strsplit(cons[i], "_")[[1]]
    erg[i] = paste(tmp[3], tmp[4])
  }
  return(erg)
}

#
# Returning species
#
get_species <- function(cons){
  erg = vector()
  for(i in 1:length(cons)){
    tmp = strsplit(cons[i], "_")[[1]]
    erg[i] = tmp[1]
  }
  return(erg)
}

#
# Returning tissues
#
get_tissues <- function(cons){
  erg = vector()
  for(i in 1:length(cons)){
    tmp = strsplit(cons[i], "_")[[1]]
    erg[i] = tmp[2]
  }
  return(erg)
}

#
# transfers a count matrix to data.frame containing log-counts
#
log.df <- function(counts, methode){
  df.x = vector()
  df.methode = vector()
  df.sample = vector()
  for(x in 1:length(counts[1,])){
    df.x = c(df.x, log(counts[,x] + 1))
    df.methode = c(df.methode, rep(methode, length(counts[,x])))
    df.sample = c(df.sample, rep(colnames(counts)[x], length(counts[,x])))
  }
  df.sample = extend_names(names = df.sample,
                           snakemake = snakemake,
                           parts = mode[[2]])
  erg = data.frame(x = df.x,
                   methode = df.methode,
                   sample = df.sample)
  return(erg)
}

#
# Getting data for plotting RLE-Plots
#
no.plotRLE <- function(data){
  file.name = paste0("RLE.", snakemake@params$mode, ".tmp.pdf")
  pdf(file.name)
  erg <- plotRLE(data, outline = F)
  dev.off()
  return(erg)
}

#
# Getting data for plotting RLE-Plots, separated by tissue/species
#
no.plotRLE.parts <- function(data){
  parts = vector()
  if(snakemake@params$mode == "species"){
    parts = get_tissues(colnames(data))
  } else {
    parts = get_species(colnames(data))
  }
  u.parts = unique(parts)
  erg = matrix(ncol = 0, nrow = length(data[,1]))
  for(x in 1:length(u.parts)){
    poses = which(parts == u.parts[x])
    if(length(poses) < 2)
      next
    subdata = data[,poses]
    file.name = paste0("RLE.", snakemake@params$mode, ".tmp.pdf")
    pdf(file.name)
    erg <- cbind(erg, plotRLE(subdata, outline = F))
    dev.off()
    file.remove(file.name)
  }
  return(erg)
}

#
# plotting ecdf
#
plotEcdf <- function(counts, output, long_spec, snakemake) {
  library(wesanderson)
  png(output, width = 720)
  layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(3, 1))
  par(mar=c(8.1,4.1,4.1,2.1))
  cols = wes_palette("Zissou1", length(counts[1,]), type = "continuous")
  plot(ecdf(log(counts[,1] + 1)), 
       main = paste("Cumulative distribution of", long_spec), col = cols[1])
  for(k in 2:length(counts[1,])){
    lines(ecdf(log(counts[,k] + 1)), col = cols[k], type = "l")
  }
  par(mar=c(5.1, 0, 4.1, 1.1), las = 2)
  plot(1, type="n", axes=F, xlab="", ylab="")
  leg = legend("topleft", 
               legend = extend_names(colnames(counts),
                                     snakemake = snakemake,
                                     parts = c(2,3)),
               pch = 22,
               col = cols,
               cex = 1, pt.bg = cols)
  dev.off()
}

#
# Methode for txi to DGEList
#
prep.DGEList <- function(txi, counts = txi$counts){
  cts <- counts
  normMat <- txi$length
  normMat <- normMat/exp(rowMeans(log(normMat)))
  o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
  groups = extract_cons(colnames(txi$counts))
  y <- DGEList(cts, group = factor(groups))
  y$offset <- t(t(log(normMat)) + o)
  y = calcNormFactors(y, method = "TMM")
  return(y)
}

set_label_levels <- function(vec){
  tmp = as.numeric(gsub(",", ".", vec, fixed = T))
  levels = tmp
  levels[which(is.na(tmp))] = "a"
  levels[which(tmp <= 100)] = "b"
  levels[which(tmp < 90)] = "c"
  levels[which(tmp < 70)] = "d"
  return(levels)
}

#
# matrix to data.frame containing geometric mean
#
sum.to.df <- function(counts, methode){
  frame.methode = vector()
  frame.counts = vector()
  for (y in 1:length(counts[,1])){
    frame.methode[length(frame.methode) + 1] = methode
    frame.counts[length(frame.counts) + 1] = prod(counts[y,])^(1/length(counts[y,]))
  }
  frame.counts = log(frame.counts + 1)
  erg <- data.frame(methode = frame.methode,
                    logcounts = frame.counts)
  return(erg)
}

#
# matrix to data.frame
#
to.df <- function(counts, methode){
  frame.sample = vector()
  frame.methode = vector()
  frame.counts = vector()
  names = colnames(counts)
  names = extend_names(names = names,
                       snakemake = snakemake,
                       parts = mode[[2]])
  for(x in 1:length(counts[1,])){
    for (y in 1:length(counts[,1])){
      frame.sample[length(frame.sample) + 1] = names[x]
      frame.methode[length(frame.methode) + 1] = methode
      frame.counts[length(frame.counts) + 1] = counts[y,x]
    }
  }
  erg <- data.frame(sample = frame.sample,
                    methode = frame.methode,
                    logcounts = frame.counts)
  return(erg)
}