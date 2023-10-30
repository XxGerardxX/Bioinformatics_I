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


peaks.AH22941 
peaks.AH22941 # contains peak information. compare to 
peaks.AH22942 
peaks.AH22942
peaks.AH22943 
peaks.AH22943
peaks.AH22944 
peaks.AH22944
peaks.AH22945 
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



#############





# Ensure that both granges_from_csv and peaks_range_granges are GRanges
granges_from_csv <- as(granges_from_csv, "GRanges")


# Now, you can use countOverlaps to check for overlaps
overlap_counts <- countOverlaps(peaks_range_granges, granges_from_csv)

# Create an overlap vector (1 for overlaps, 0 for no overlaps)
overlap_vector <- ifelse(overlap_counts > 0, 1, 0)

overlap_vector


peaks.AH22941
granges_from_csv


#ranges object with start end and with
peaks_range <- ranges(peaks.AH22941)
peaks_range


findOverlaps(granges_from_csv,peaks.AH22941)


str(peaks.AH22941)
# Convert the IRanges object to a GRanges object with single positions




length(encodeFiles.HepG2) # file length is 48
seqnames(encodeFiles.HepG2)


# Create two sequences
seq1 <- 22941:22975
seq2 <- 23305:23317




write.csv(data.frame(start = start(peaks_range), end = end(peaks_range), width = width(peaks_range)), "peaks_data.csv", row.names = FALSE)

