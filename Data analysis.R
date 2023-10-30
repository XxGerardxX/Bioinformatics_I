library(SummarizedExperiment)
library(writexl)
lihc <- readRDS("LIHC_methyl450k_subset_summexp.rds")
lihc
rowData(lihc)







chromosome <- rowData(lihc)$Chromosome 
genomic_coord <- rowData(lihc)$Genomic_Coordinate
assay(lihc)

rowd <- rowData(lihc)
assay <- assay(lihc)


write.csv(assay, "assay.csv")
rowd
assay
