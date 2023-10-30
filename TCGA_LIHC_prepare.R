# R4.2.1@ux-p-bicalc
#
# See /data/home/pdmoerland/moerland/TCGA/Bioinformatics-I

library(curatedTCGAData)
library(SummarizedExperiment)
library(TCGAutils)

curatedTCGAData(diseaseCode = "*", assays = "*", dry.run = TRUE, version = "2.0.1")

# LIHC
curatedTCGAData(diseaseCode = "LIHC", assays = "methyl*", dry.run = TRUE, version = "2.0.1")
lihc <- curatedTCGAData("LIHC", "methyl*", FALSE, version = "2.0.1")
saveRDS(lihc,file="LIHC_methyl_curatedTCGAData.rds")

lihc <- readRDS("LIHC_methyl_curatedTCGAData.rds")

# What sample types are present in the data?
sampleTables(lihc)
# Select 'Primary Solid Tumor' samples
lihc.tm <- TCGAsplitAssays(lihc, "01")
# Sample only 50 patients
set.seed(11)
lihc.subset <- lihc.tm[,sample(1:ncol(lihc.tm[[1]]),50)]
saveRDS(lihc.subset,"LIHC_methyl450k_subset.rds")

lihc <- readRDS("LIHC_methyl450k_subset.rds")

# Switch from a MultiAssayExperiment to a SummarizedExperiment
lihc.summexp <- lihc[[1]]

# Remove all rows that contain at least one missing value
lihc.summexp <- lihc.summexp[rowSums(is.na(assay(lihc.summexp)))==0,]
# Remove features that are not linked to a chromosome
lihc.summexp <- lihc.summexp[!is.na(rowData(lihc.summexp)$Chromosome),]

table((rowData(lihc.summexp)$Chromosome))
saveRDS(lihc.summexp,"LIHC_methyl450k_subset_summexp.rds")

# Here it starts
lihc.summexp <- readRDS("LIHC_methyl450k_subset_summexp.rds")

# They switched to DelayedMatrix and this only works if you also have access to the cached files
assay(lihc.summexp)

# So let's convert it to a numeric matrix
lihc.summexp <- SummarizedExperiment(assays=list(counts=as.matrix(assay(lihc.summexp))), rowData=rowData(lihc.summexp), colData=colData(lihc.summexp))
saveRDS(lihc.summexp,"LIHC_methyl450k_subset_summexp.rds")

# Here it starts again
lihc <- readRDS("LIHC_methyl450k_subset_summexp.rds")

