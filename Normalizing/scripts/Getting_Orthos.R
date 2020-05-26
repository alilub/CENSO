#
# Function for downloading othologe genes form ensemble
#

library(refGenome)
library("biomaRt")

gtf <- ensemblGenome()
link <- paste0("../", "Human", "/References/GTF/")
ref <- list.files(link, full.names = T)
read.gtf(gtf, ref)

# Herauslesen der gene_id

tmp = getGeneTable(gtf)
gene_ids = unique(tmp$gene_id)

ensembl=useMart("ensembl")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl", host = "http://apr2018.archive.ensembl.org")

orthos = getBM(c("ensembl_gene_id", 
                 "ppaniscus_homolog_orthology_type", "ppaniscus_homolog_ensembl_gene",
                 "ptroglodytes_homolog_orthology_type", "ptroglodytes_homolog_ensembl_gene",
                 "ggallus_homolog_orthology_type", "ggallus_homolog_ensembl_gene",
                 "mmusculus_homolog_orthology_type", "mmusculus_homolog_ensembl_gene",
                 "ggorilla_homolog_orthology_type", "ggorilla_homolog_ensembl_gene"
                 ),
               filters = c("ensembl_gene_id", 
                           "with_ppaniscus_homolog", "with_ptroglodytes_homolog", "with_ggallus_homolog",
                           "with_mmusculus_homolog", "with_ggorilla_homolog", "with_mmulatta_homolog",
                           "with_mdomestica_homolog", "with_pabelii_homolog", "with_oanatinus_homolog"),
               values = list(gene_ids, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
               mart = ensembl)

tmp = intersect(which(orthos[,2] == "ortholog_one2one"),
                which(orthos[,4] == "ortholog_one2one"))
tmp = intersect(tmp, which(orthos[,6] == "ortholog_one2one"))
tmp = intersect(tmp, which(orthos[,8] == "ortholog_one2one"))
tmp = intersect(tmp, which(orthos[,10] == "ortholog_one2one"))
orthos = orthos[tmp,]

save.image("Orthos_tmp.RData")

orthos_2 = getBM(c("ensembl_gene_id",
                   "mdomestica_homolog_orthology_type", "mdomestica_homolog_ensembl_gene",
                   "pabelii_homolog_orthology_type", "pabelii_homolog_ensembl_gene",
                   "oanatinus_homolog_orthology_type", "oanatinus_homolog_ensembl_gene",
                   "mmulatta_homolog_orthology_type", "mmulatta_homolog_ensembl_gene"
                   ),
                 filters = c("ensembl_gene_id", 
                             "with_ppaniscus_homolog", "with_ptroglodytes_homolog", "with_ggallus_homolog",
                             "with_mmusculus_homolog", "with_ggorilla_homolog", "with_mmulatta_homolog",
                             "with_mdomestica_homolog", "with_pabelii_homolog", "with_oanatinus_homolog"),
                 values = list(gene_ids, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                 mart = ensembl)

tmp = intersect(which(orthos_2[,2] == "ortholog_one2one"),
                which(orthos_2[,4] == "ortholog_one2one"))
tmp = intersect(tmp, which(orthos_2[,6] == "ortholog_one2one"))
tmp = intersect(tmp, which(orthos_2[,8] == "ortholog_one2one"))
orthos_2 = orthos_2[tmp,]

genes = intersect(orthos[,1], orthos_2[,1])
poses = match(genes, orthos[,1])
orthos = orthos[poses,]

poses_2 = match(genes, orthos_2[,1])
orthos_2 = orthos_2[poses_2,]

all_orthos = cbind(orthos, orthos_2)
write.csv(all_orthos, sep = ";", quote = FALSE, file = "orthos.csv")
