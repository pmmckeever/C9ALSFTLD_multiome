library(Sierra)
library(dplyr)
library(stringr)
workingdir <- '/data1/APA/Paul_ALS_Data/sierra_out_2/'
setwd(workingdir)
reference.file <- '/home/aiden/data/refgenome/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
root <- '/data1/APA/Paul_ALS_Data/bams_in/'
samples <- read.table(paste0(root,'metadata_CTpaul_ALSpaul.txt'), header = T)
peak.output.file <- c(paste0(samples$name, '_FindPeak_out'))
L <- dim(samples)[1]
print('Step 1: Samples are read and now performing the peak calling per sample')
#for (row in 1:L){
#    FindPeaks(output.file = peak.output.file[row],      # output filename
#    gtf.file = reference.file,                   # gene model as a GTF file
#    bamfile = paste0(root,samples[row,]$bam),                   # BAM alignment filename.
#    junctions.file = paste0(root,samples[row,]$junctions),     # BED filename of splice junctions exising in BAM file.
#    ncores = 32)
#}

print('Step 2: peak calling per sample is done now mergeing the peaks')
#peak.dataset.table = data.frame(Peak_file = peak.output.file,
#  Identifier = samples$name,
#  stringsAsFactors = FALSE)
peak.merge.output.file = "ALS_paul_merged_peaks.txt"
#MergePeakCoordinates(peak.dataset.table, output.file = peak.merge.output.file, ncores = 32)
####
print('Step 3: merging is done now counting the per peak, this can take a while')

count.dir <- paste0(samples$name, "_sierra_counts")
for (row in 1:L){
    print(count.dir[row])
    CountPeaks(peak.sites.file = peak.merge.output.file,
    gtf.file = reference.file,
    bamfile = paste0(root,samples[row,]$bam),
    whitelist.file = paste0(root,samples[row,]$barcodes),
    output.dir = count.dir[row],
    countUMI = TRUE,
    ncores = 32)

}

print('Step 4: counting is done now aggeragating the peaks for all samples')
count.dir <- paste0(samples$name, "_sierra_counts")
barcode_extensions <- data.frame(samples$name)
colnames(barcode_extensions) <- c('ext')
barcode_extensions <- paste0("_",barcode_extensions$ext)

out.dir <- "ALS_paul_Sierra_aggregate"
#Now aggregate the counts for all datasets
AggregatePeakCounts(peak.sites.file = peak.merge.output.file,
                    count.dirs = count.dir,
                    exp.labels = barcode_extensions,
                    output.dir = out.dir)

genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
print('Aggeration is done; now annotating the peaks')
AnnotatePeaksFromGTF(peak.sites.file = peak.merge.output.file,
                     gtf.file = reference.file,
                     output.file = "ALS_paul_peaks_annotated.txt",
                     genome = genome)
#print('Last step 5: reading and writing the annotated peak counts')
#peak.counts <- ReadPeakCounts(data.dir = 'Sierra_out_1/Kapmann_Sierra_aggregate/')
#Read in peak annotations
#peak.annotations <- read.table("Kapmann_peaks_annotated.txt",
#                               header = TRUE,
#                               sep = "\t",
#                               row.names = 1,
#                               stringsAsFactors = FALSE)
#write.table(peak.annotations, file=paste0(workingdir, 'peak.annotations.tsv'), sep='\t')
