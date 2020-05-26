#
# Plotting the scale factors
#
source("scripts/Functions.R")
suppressPackageStartupMessages(library("viridis"))
suppressPackageStartupMessages(library(ggplot2))
#
# for species
#
methods = as.vector(read.csv(snakemake@config$methods, sep = "\t")[,1])
methods = c("raw", methods)
save.image(paste0("images/", snakemake@rule, ".RData"))

#
# transfering the input
#
in.DESeq= snakemake@input$all_deseq
in.TMM = snakemake@input$all_tmm
in.RPKM_TMM = snakemake@input$all_rpkm_tmm
in.RPKM_DESeq = snakemake@input$all_rpkm_deseq
in.TPM_TMM = snakemake@input$all_tpm_tmm
in.TPM_DESeq = snakemake@input$all_tpm_deseq
output = snakemake@output$all
ps.names = c(1,3)
ps.part = c(2)

#
# Plotting
#
for (i in 1:length(in.DESeq)) {
  # Data DESeq
  load(in.DESeq[i])
  sf.deseq = sf
  names(sf.deseq) = NULL
  # Data TMM
  load(in.TMM[i])
  sf.tmm = sf
  sf.names = extend_names(names = names(sf.tmm), 
                          snakemake =  snakemake, 
                          parts = ps.names)
  sf.part = extend_names(names = names(sf.tmm), 
                         snakemake =  snakemake, 
                         parts = ps.part)
  sf.part = rep(sf.part, 6)
  names(sf.tmm) = NULL
  # Data RPKM_DESeq
  load(in.RPKM_DESeq[i])
  sf.rpkm_deseq = sf
  names(sf.rpkm_deseq) = NULL
  # Data RPKM_TMM
  load(in.RPKM_TMM[i])
  sf.rpkm_tmm = sf
  names(sf.rpkm_tmm) = NULL
  # Data TPM_DESeq
  load(in.TPM_DESeq[i])
  sf.tpm_deseq = sf
  names(sf.tpm_deseq) = NULL
  # Data TPM_TMM
  load(in.TPM_TMM[i])
  sf.tpm_tmm = sf
  names(sf.tpm_tmm) = NULL
  #
  # to data.frame
  #
  sf.values = c(sf.tmm, 
                sf.deseq,
                sf.rpkm_tmm,
                sf.rpkm_deseq,
                sf.tpm_tmm,
                sf.tpm_deseq)
  sf.methode = c(rep("TMM", length(sf.tmm)),
                 rep("DESeq", length(sf.deseq)),
                 rep("RPKM-TMM", length(sf.rpkm_tmm)),
                 rep("RPKM-DESeq", length(sf.rpkm_deseq)),
                 rep("TPM-TMM", length(sf.tpm_tmm)),
                 rep("TPM-DESeq", length(sf.tpm_deseq)))
  sf.methode = factor(sf.methode, levels = unique(sf.methode))
  sf.sample = rep(sf.names, 6)
  sf = data.frame(value = sf.values,
                  methode = sf.methode,
                  sample = sf.sample)
  for(y in 1:length(unique(sf.part))){
    part = sort(unique(sf.part))[y]
    pos = which(sf.part == part)
    data = sf[pos,]
    bp <- ggplot(data=data, aes(x=sample, y=value + 1, fill=sample)) +
      theme_linedraw()+
      geom_bar(stat="identity", position=position_dodge(), color = "black")+
      geom_text(aes(label = round(value,2)), size = 2.5, vjust=-0.2, position=position_dodge(1))+
      scale_fill_manual(values = get_colour_met(length(pos)/6), name = "sample")+
      theme(axis.text.x=element_text(angle = -90, hjust = 0))+
      ggtitle(paste("Scale factors of", part))+
      scale_y_continuous(name = "log10 scale of value + 1",
                         trans = "log10")+
      facet_grid(methode~.)
    bp + geom_hline(yintercept = 3, color='red', linetype = "dashed")
    png(output[y],
        width = 1440,
        height = 720,
        res = 120)
    par(mar=c(9.1,4.1,4.1,2.1))
    print(bp)
    dev.off()
  }
} 