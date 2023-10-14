library(SummarizedExperiment)
library(BiocManager)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
?IlluminaHumanMethylation450kanno.ilmn12.hg19

illumina_data <- data(IlluminaHumanMethylation450kanno.ilmn12.hg19)


table(Islands.UCSC[,2])

# island data from shelf=, shore or island make into 3 groups
write.csv(Islands.UCSC, file = "island_data.csv", row.names = TRUE)


# print(unique(Islands.UCSC$Relation_to_Island))

# In the refgene_group is the gene body, body, promotor and intergenic region
print(Other)
Other_data <- (Other$UCSC_RefGene_Group)

write.csv(Other, file = "Other.csv", row.names = TRUE)

#snps illumina

SNP_ill <- data("SNPs.Illumina")
write.csv(SNPs.Illumina, file = "SNPs.csv", row.names = TRUE)



