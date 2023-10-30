if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")



BiocManager::install("AnnotationHub")
BiocManager::install("GenomicRanges")

library(AnnotationHub)
library(GenomicRanges)
library(data.table)
library(IRanges)

ah <- AnnotationHub()
encodeFiles.HepG2 <- query(ah,c("ENCODE","HepG2","TfbsSydh|BroadHistone"))
encodeFiles.HepG2
# AnnotationHub with 48 records
# snapshotDate(): 2020-04-27
# $dataprovider: UCSC
# $species: Homo sapiens
# $rdataclass: GRanges
# additional mcols(): taxonomyid, genome, description, coordinate_1_based, maintainer, rdatadateadded, preparerclass, tags,
#   rdatapath, sourceurl, sourcetype 
# retrieve records with, e.g., 'object[["AH22941"]]' 
genomic_distance_data <- read.csv('brca.csv')



peaks.AH22941 <- ah[["AH22941"]]
peaks.AH22941 # contains peak information. compare to 
peaks.AH22942 <-ah[["AH22942"]]
peaks.AH22942
peaks.AH22943 <-ah[["AH22943"]]
peaks.AH22943
peaks.AH22944 <-ah[["AH22944"]]
peaks.AH22944
peaks.AH22945 <-ah[["AH22945"]]
peaks.AH22946
peaks.AH22947
peaks.AH22948
peaks.AH22949
peaks.AH22950
peaks.AH22951
peaks.AH22952
peaks.AH22953
peaks.AH22954
peaks.AH22955
peaks.AH22956
peaks.AH22957
peaks.AH22958
peaks.AH22959
peaks.AH22960
peaks.AH22961
peaks.AH22962
peaks.AH22963
peaks.AH22964
peaks.AH22965
peaks.AH22966
peaks.AH22967
peaks.AH22968
peaks.AH22969
peaks.AH22970
peaks.AH22971
peaks.AH22972
peaks.AH22973
peaks.AH22974
peaks.AH22975

# Access the peaks objects from AH23305 to AH23317
peaks.AH23305 
peaks.AH23306 
peaks.AH23307 
peaks.AH23308 
peaks.AH23309 
peaks.AH23310 
peaks.AH23311 
peaks.AH23312 
peaks.AH23313 
peaks.AH23314 
peaks.AH23315 
peaks.AH23316 
peaks.AH23317 






############Loops over the data to create peak data#############


# Create a list of AH object names you want to access
ah_object_names <- sprintf("AH%05d", 22941:22975)

# Create objects with the names peaks.AHXXXX
for (ah_name in ah_object_names) {
  assign(paste("peaks.", ah_name, sep = ""), ah[[ah_name]])
}


# Create a list of AH object names you want to access
ah_object_names_2 <- sprintf("AH%05d", 23305:23317)

# Create objects with the names peaks.AHXXXX
for (ah_name in ah_object_names_2) {
  assign(paste("peaks.", ah_name, sep = ""), ah[[ah_name]])
}


################################################################

# Step 1: Read data from "brca.csv" and create a GRanges object
# Replace 'brca.csv' with the actual path to your CSV file
brca_data <- brca
# Assuming 'brca_data' is your data.table containing genomic coordinate information
# Create a GRanges object
granges_from_data_table <- GRanges(
  seqnames = sub("^", "chr", brca_data$Chromosome)
,  # Replace 'seqnames' with the appropriate column name
  ranges = IRanges(
    start = as.integer(brca_data$Genomic_Coordinate
),  # Replace 'start' with the appropriate column name
    end = as.integer(brca_data$Genomic_Coordinate
)       # Replace 'end' with the appropriate column name
  )
)


# Convert granges_from_csv to a GRanges object with single positions

granges_from_data_table


# Step 2: Loop through each peaks object and compare overlaps
ah_peak_files <- c(
  "peaks.AH22941", "peaks.AH22942", "peaks.AH22943", "peaks.AH22944", "peaks.AH22945",
  "peaks.AH22946", "peaks.AH22947", "peaks.AH22948", "peaks.AH22949", "peaks.AH22950",
  "peaks.AH22951", "peaks.AH22952", "peaks.AH22953", "peaks.AH22954", "peaks.AH22955",
  "peaks.AH22956", "peaks.AH22957", "peaks.AH22958", "peaks.AH22959", "peaks.AH22960",
  "peaks.AH22961", "peaks.AH22962", "peaks.AH22963", "peaks.AH22964", "peaks.AH22965",
  "peaks.AH22966", "peaks.AH22967", "peaks.AH22968", "peaks.AH22969", "peaks.AH22970",
  "peaks.AH22971", "peaks.AH22972", "peaks.AH22973", "peaks.AH22974", "peaks.AH22975",
  "peaks.AH23305", "peaks.AH23306", "peaks.AH23307", "peaks.AH23308", "peaks.AH23309",
  "peaks.AH23310", "peaks.AH23311", "peaks.AH23312", "peaks.AH23313", "peaks.AH23314",
  "peaks.AH23315", "peaks.AH23316", "peaks.AH23317"
)





# Convert the IRanges object to a GRanges object with single positions
peaks <- matrix(0, nrow(brca_data), length(encodeFiles.HepG2))
colnames(peaks) <- names(encodeFiles.HepG2)
for (i in names(encodeFiles.HepG2)){
  peaks[, i] <- countOverlaps(granges_from_data_table, ah[[i]])
}
peaks[peaks>1] <- 1
colnames(peaks) <- sub(".narrowPeak.gz|.broadPeak.gz","", sub("wgEncode(\\w*)Hepg2","",encodeFiles.HepG2$title))



write.csv(peaks, file = "48_peaks_allignment.csv", row.names = TRUE)










