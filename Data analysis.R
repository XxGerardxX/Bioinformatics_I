
library(writexl)
chromosome <- rowData(brca)$Chromosome 
genomic_coord <- rowData(brca)$Genomic_Coordinate
assay(brca)

rowd <- rowData(brca)
assay <- assay(brca)


write.csv(assay, "assay.csv")
rowd
assay
