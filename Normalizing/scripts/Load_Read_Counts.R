#
# Script for transfering salmon output into txi objects (one per species)
#
# Carful, can only handle two technical replicates!
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(edgeR))
#
# Extracting the species
#
save.image(paste0("images/", snakemake@rule, ".RData"))

spec <- vector()
for(i in 1:length(snakemake@output))
{
  tmp = strsplit(snakemake@output[[i]], "/")[[1]][4]
  spec[i]= strsplit(tmp, ".", fixed = T)[[1]][1]
}
#
# Creating a txi object per species
#
for(s in 1:length(spec)){
    species = spec[s]
    print(species)
      link <- paste0("../../../hiwi/Brawand/data/", species, "/References/GTF/")
      ref <- list.files(link, full.names = T)
      print(ref)
      ens <- makeTxDbFromGFF(ref)
      k <- keys(ens, keytype = "TXNAME")
      tx2gene <- select(ens, k, "GENEID", "TXNAME")
    #
    # Getting and naming the files
    #
    tmp <- unlist(snakemake@input)
    n <- vector()
    files <- vector()
    #
    # filtering the files
    #
    for(i in 1:length(tmp)){
      if(length(grep(species, tmp[i])) > 0){
        files[length(files) + 1] = tmp[i]
      }
    }
    for(i in 1:length(files))
      n[i] <- strsplit(files[i], "/")[[1]][12]
    names(files) <- n
    files  = sub(pattern = "/scratch/user/acfhm/", replacement = "../../../", x = files, fixed = T)
    #
    # Import read counts
    #
    txi <- tximport(files, type = "salmon", tx2gene = tx2gene, dropInfReps = T)
    names(txi)
    #
    # sum up technical replicates
    #
    cols = colnames(txi$counts)
    sam = cols
    for(j in 1:length(sam)){
      tmp = strsplit(sam[j], "_")[[1]]
      sam[j] = paste(tmp[1], tmp[2], tmp[3], tmp[4], sep = "_")
    }
    txi$counts = sumTechReps(txi$counts, ID = sam)
    txi$abundance = sumTechReps(txi$abundance, ID = sam)
    techRep = unique(sam[duplicated(sam)])
    techRep.pos = match(techRep, colnames(txi$abundance))
    txi$abundance[,techRep.pos] = txi$abundance[,techRep.pos]/2
    txi$length = txi$length[,match(unique(sam), sam)]
    save(txi, file = paste0("data/counts/raw/", species, ".RData"))
}







