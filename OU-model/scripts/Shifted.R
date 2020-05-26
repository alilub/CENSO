methods = as.vector(read.csv(snakemake@config$methods_all, sep = "\t")[,1])

sig.genes = list()
for(i in 1:length(snakemake@input)){
  theta.LRT = as.vector(t(read.table(snakemake@input[[i]])))
  gene.ids = as.vector(read.table(
    paste0("data/eve/all/", methods[i], "/", snakemake@wildcards$tis, "/",snakemake@wildcards$tis,".dat"), skip = 1)$V1)
  # Calculate p-values
  ptheta = pchisq(theta.LRT, 1)
  # Correct p-values
  ptheta.c = p.adjust(ptheta, method = "BY")
  sig = which(ptheta.c < 0.05)
  gene.ids.sig = gene.ids[sig]
  sig.genes[[i]] = gene.ids.sig
}

erg = matrix(nrow = length(methods),
             ncol = length(methods), 
             dimnames = list(methods, methods))

erg.rel = erg = matrix(nrow = length(methods),
                       ncol = length(methods), 
                       dimnames = list(methods, methods))
for(i in 1:length(methods)){
  for(j in 1:length(methods)){
    erg[i,j] = length(intersect(sig.genes[[i]], sig.genes[[j]]))
    erg.rel[i,j] = erg[i,j] / length(unique(c(sig.genes[[i]], sig.genes[[j]])))
  }
}

png(filename =  snakemake@output[[1]], 
    width = 960,
    height = 960)
pheatmap::pheatmap(erg, display_numbers = T, cex = 1.5, 
                   cluster_cols = FALSE,
                   cluster_rows = FALSE)
graphics.off()

png(filename =  snakemake@output[[2]], 
    width = 960,
    height = 960)
pheatmap::pheatmap(erg.rel, display_numbers = T, cex = 1.5, 
                   cluster_cols = FALSE,
                   cluster_rows = FALSE)
graphics.off()