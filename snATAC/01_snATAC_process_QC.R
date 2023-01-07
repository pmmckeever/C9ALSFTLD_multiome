library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(future)
set.seed(1234)

#control_cases
seur01 <- read.table(
  file = "ALS_snATAC/CTRL1_snATAC/peaks.bed",
  col.names = c("chr", "start", "end")
)
seur02 <- read.table(
  file = "ALS_snATAC/CTRL2_snATAC/peaks.bed",
  col.names = c("chr", "start", "end")
)
seur03 <- read.table(
  file = "ALS_snATAC/CTRL3_snATAC/peaks.bed",
  col.names = c("chr", "start", "end")
)
seur04 <- read.table(
  file = "ALS_snATAC/CTRL4_snATAC/peaks.bed",
  col.names = c("chr", "start", "end")
)


#C9_noALS_noFTLD
seur70 <- read.table(
  file = "ALS_snATAC/C9noALSnoFTLD_snATAC/peaks.bed",
  col.names = c("chr", "start", "end")
)

#C9ALSFTLD
seur11 <- read.table(
  file = "ALS_snATAC/C9ALSFTLD1_snATAC/peaks.bed",
  col.names = c("chr", "start", "end")
)
seur12 <- read.table(
  file = "ALS_snATAC/C9ALSFTLD2_snATAC/peaks.bed",
  col.names = c("chr", "start", "end")
)
seur13 <- read.table(
  file = "ALS_snATAC/C9ALSFTLD3_snATAC/peaks.bed",
  col.names = c("chr", "start", "end")
)
seur14 <- read.table(
  file = "ALS_snATAC/C9ALSFTLD5_snATAC/peaks.bed",
  col.names = c("chr", "start", "end")
)
seur15 <- read.table(
  file = "ALS_snATAC/C9ALSFTLD6_snATAC/peaks.bed",
  col.names = c("chr", "start", "end")
)
seur16 <- read.table(
  file = "ALS_snATAC/C9ALSFTLD7_snATAC/peaks.bed",
  col.names = c("chr", "start", "end")
)

#C9ALSnoFTLD
seur21 <- read.table(
  file = "ALS_snATAC/C9ALSnoFTLD1_snATAC/peaks.bed",
  col.names = c("chr", "start", "end")
)
seur22 <- read.table(
  file = "ALS_snATAC/C9ALSnoFTLD2_snATAC/peaks.bed",
  col.names = c("chr", "start", "end")
)
seur23 <- read.table(
  file = "ALS_snATAC/C9ALSnoFTLD3_snATAC/peaks.bed",
  col.names = c("chr", "start", "end")
)

#sALSnoFTLD
seur1 <- read.table(
  file = "ALS_snATAC/sALSnoFTLD1_snATAC/peaks.bed",
  col.names = c("chr", "start", "end")
)
seur2 <- read.table(
  file = "ALS_snATAC/sALSnoFTLD2_snATAC/peaks.bed",
  col.names = c("chr", "start", "end")
)
seur3 <- read.table(
  file = "ALS_snATAC/sALSnoFTLD3_snATAC/peaks.bed",
  col.names = c("chr", "start", "end")
)
seur4 <- read.table(
  file = "ALS_snATAC/sALSnoFTLD4_snATAC/peaks.bed",
  col.names = c("chr", "start", "end")
)
seur5 <- read.table(
  file = "ALS_snATAC/sALSnoFTLD5_snATAC/peaks.bed",
  col.names = c("chr", "start", "end")
)
seur6 <- read.table(
  file = "ALS_snATAC/sALSnoFTLD6_snATAC/peaks.bed",
  col.names = c("chr", "start", "end")
)
seur7 <- read.table(
  file = "ALS_snATAC/sALSnoFTLD7_snATAC/peaks.bed",
  col.names = c("chr", "start", "end")
)
seur8 <- read.table(
  file = "ALS_snATAC/sALSnoFTLD8_snATAC/peaks.bed",
  col.names = c("chr", "start", "end")
)

# convert to genomic ranges
gr.01 <- makeGRangesFromDataFrame(seur01)
gr.02 <- makeGRangesFromDataFrame(seur02)
gr.03 <- makeGRangesFromDataFrame(seur03)
gr.04 <- makeGRangesFromDataFrame(seur04)
gr.70 <- makeGRangesFromDataFrame(seur70)
gr.11 <- makeGRangesFromDataFrame(seur11)
gr.12 <- makeGRangesFromDataFrame(seur12)
gr.13 <- makeGRangesFromDataFrame(seur13)
gr.14 <- makeGRangesFromDataFrame(seur14)
gr.15 <- makeGRangesFromDataFrame(seur15)
gr.16 <- makeGRangesFromDataFrame(seur16)
gr.21 <- makeGRangesFromDataFrame(seur21)
gr.22 <- makeGRangesFromDataFrame(seur22)
gr.23 <- makeGRangesFromDataFrame(seur23)
gr.1 <- makeGRangesFromDataFrame(seur1)
gr.2 <- makeGRangesFromDataFrame(seur2)
gr.3 <- makeGRangesFromDataFrame(seur3)
gr.4 <- makeGRangesFromDataFrame(seur4)
gr.5 <- makeGRangesFromDataFrame(seur5)
gr.6 <- makeGRangesFromDataFrame(seur6)
gr.7 <- makeGRangesFromDataFrame(seur7)
gr.8 <- makeGRangesFromDataFrame(seur8)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = c(gr.01, gr.02, gr.03, gr.04, gr.70, gr.11, gr.12, gr.13, gr.15, gr.16, gr.17, gr.21, gr.22, gr.23,
                               gr.1, gr.2, gr.3, gr.4, gr.5, gr.6, gr.7, gr.8))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

# load metadata
md.01 <- read.table(
  file = "ALS_snATAC/CTRL1_snATAC/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.02 <- read.table(
  file = "ALS_snATAC/CTRL2_snATAC/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.03 <- read.table(
  file = "ALS_snATAC/CTRL3_snATAC/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.04 <- read.table(
  file = "ALS_snATAC/CTRL4_snATAC/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.70 <- read.table(
  file = "ALS_snATAC/C9noALSnoFTLD_snATAC/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]


md.11 <- read.table(
  file = "ALS_snATAC/C9ALSFTLD1_snATAC/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.12 <- read.table(
  file = "ALS_snATAC/C9ALSFTLD2_snATAC/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.13 <- read.table(
  file = "ALS_snATAC/C9ALSFTLD3_snATAC/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.14 <- read.table(
  file = "ALS_snATAC/C9ALSFTLD5_snATAC/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.15 <- read.table(
  file = "ALS_snATAC/C9ALSFTLD6_snATAC/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.16 <- read.table(
  file = "ALS_snATAC/C9ALSFTLD7_snATAC/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.21 <- read.table(
  file = "ALS_snATAC/C9ALSnoFTLD1_snATAC/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.22 <- read.table(
  file = "ALS_snATAC/C9ALSnoFTLD2_snATAC/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.23 <- read.table(
  file = "ALS_snATAC/C9ALSnoFTLD3_snATAC/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.1 <- read.table(
  file = "ALS_snATAC/sALSnoFTLD1_snATAC/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.2 <- read.table(
  file = "ALS_snATAC/sALSnoFTLD2_snATAC/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.3 <- read.table(
  file = "ALS_snATAC/sALSnoFTLD3_snATAC/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.4 <- read.table(
  file = "ALS_snATAC/sALSnoFTLD4_snATAC/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.5 <- read.table(
  file = "ALS_snATAC/sALSnoFTLD5_snATAC/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.6 <- read.table(
  file = "ALS_snATAC/sALSnoFTLD6_snATAC/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.7 <- read.table(
  file = "ALS_snATAC/sALSnoFTLD7_snATAC/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

md.8 <- read.table(
  file = "ALS_snATAC/sALSnoFTLD8_snATAC/singlecell.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ]

# perform an initial filtering of low count cells
md.01 <- md.01[md.01$passed_filters > 500, ]
md.02 <- md.02[md.02$passed_filters > 500, ]
md.03 <- md.03[md.03$passed_filters > 500, ]
md.04 <- md.04[md.04$passed_filters > 500, ]
md.70 <- md.70[md.70$passed_filters > 500, ]
md.11 <- md.11[md.11$passed_filters > 500, ]
md.12 <- md.12[md.12$passed_filters > 500, ]
md.13 <- md.13[md.13$passed_filters > 500, ]
md.14 <- md.14[md.14$passed_filters > 500, ]
md.15 <- md.15[md.15$passed_filters > 500, ]
md.16 <- md.16[md.16$passed_filters > 500, ]
md.21 <- md.21[md.21$passed_filters > 500, ]
md.22 <- md.22[md.22$passed_filters > 500, ]
md.23 <- md.23[md.23$passed_filters > 500, ]
md.1 <- md.1[md.1$passed_filters > 500, ]
md.2 <- md.2[md.2$passed_filters > 500, ]
md.3 <- md.3[md.3$passed_filters > 500, ]
md.4 <- md.4[md.4$passed_filters > 500, ]
md.5 <- md.5[md.5$passed_filters > 500, ]
md.6 <- md.6[md.6$passed_filters > 500, ]
md.7 <- md.7[md.7$passed_filters > 500, ]
md.8 <- md.8[md.8$passed_filters > 500, ]

# create fragment objects
frags.01 <- CreateFragmentObject(
  path = "ALS_snATAC/fragments/fragments.tsv.gz",
  cells = rownames(md.01)
)
frags.02 <- CreateFragmentObject(
  path = "ALS_snATAC/fragments/fragments.tsv.gz",
  cells = rownames(md.02)
)
  path = "ALS_snATAC/fragments/fragments.tsv.gz",
frags.03 <- CreateFragmentObject(
  cells = rownames(md.03)
)
frags.04 <- CreateFragmentObject(
  path = "ALS_snATAC/fragments/fragments.tsv.gz",
  cells = rownames(md.04)
)

frags.70 <- CreateFragmentObject(
  path = "ALS_snATAC/fragments/fragments.tsv.gz",
  cells = rownames(md.70)
)


frags.11 <- CreateFragmentObject(
  path = "ALS_snATAC/fragments/fragments.tsv.gz",
  cells = rownames(md.11)
)

frags.12 <- CreateFragmentObject(
  path = "ALS_snATAC/fragments/fragments.tsv.gz",
  cells = rownames(md.12)
)

frags.13 <- CreateFragmentObject(
  path = "ALS_snATAC/fragments/fragments.tsv.gz",
  cells = rownames(md.13)
)

frags.14 <- CreateFragmentObject(
  path = "ALS_snATAC/fragments/fragments.tsv.gz",
  cells = rownames(md.14)
)

frags.15 <- CreateFragmentObject(
  path = "ALS_snATAC/fragments/fragments.tsv.gz",
  cells = rownames(md.15)
)

frags.16 <- CreateFragmentObject(
  path = "ALS_snATAC/fragments/fragments.tsv.gz",
  cells = rownames(md.16)
)

frags.21 <- CreateFragmentObject(
  path = "ALS_snATAC/fragments/fragments.tsv.gz",
  cells = rownames(md.21)
)

frags.22 <- CreateFragmentObject(
  path = "ALS_snATAC/fragments/fragments.tsv.gz",
  cells = rownames(md.22)
)

frags.23 <- CreateFragmentObject(
  path = "ALS_snATAC/fragments/fragments.tsv.gz",
  cells = rownames(md.23)
)

frags.1 <- CreateFragmentObject(
  path = "ALS_snATAC/fragments/fragments.tsv.gz",
  cells = rownames(md.1)
)

frags.2 <- CreateFragmentObject(
  path = "ALS_snATAC/fragments/fragments.tsv.gz",
  cells = rownames(md.2)
)

frags.3 <- CreateFragmentObject(
  path = "ALS_snATAC/fragments/fragments.tsv.gz",
  cells = rownames(md.3)
)

frags.4 <- CreateFragmentObject(
  path = "ALS_snATAC/fragments/fragments.tsv.gz",
  cells = rownames(md.4)
)

frags.5 <- CreateFragmentObject(
  path = "ALS_snATAC/fragments/fragments.tsv.gz",
  cells = rownames(md.5)
)

frags.6 <- CreateFragmentObject(
  path = "ALS_snATAC/fragments/fragments.tsv.gz",
  cells = rownames(md.6)
)

frags.7 <- CreateFragmentObject(
  path = "ALS_snATAC/fragments/fragments.tsv.gz",
  cells = rownames(md.7)
)

frags.8 <- CreateFragmentObject(
  path = "ALS_snATAC/fragments/fragments.tsv.gz",
  cells = rownames(md.8)
)

# create count matrices
seur01.counts <- FeatureMatrix(
  fragments = frags.01,
  features = combined.peaks,
  cells = rownames(md.01)
)

seur02.counts <- FeatureMatrix(
  fragments = frags.02,
  features = combined.peaks,
  cells = rownames(md.02)
)

seur03.counts <- FeatureMatrix(
  fragments = frags.03,
  features = combined.peaks,
  cells = rownames(md.03)
)

seur04.counts <- FeatureMatrix(
  fragments = frags.04,
  features = combined.peaks,
  cells = rownames(md.04)
)

seur70.counts <- FeatureMatrix(
  fragments = frags.70,
  features = combined.peaks,
  cells = rownames(md.70)
)


seur11.counts <- FeatureMatrix(
  fragments = frags.11,
  features = combined.peaks,
  cells = rownames(md.11)
)

seur12.counts <- FeatureMatrix(
  fragments = frags.12,
  features = combined.peaks,
  cells = rownames(md.12)
)

seur13.counts <- FeatureMatrix(
  fragments = frags.13,
  features = combined.peaks,
  cells = rownames(md.13)
)

seur14.counts <- FeatureMatrix(
  fragments = frags.14,
  features = combined.peaks,
  cells = rownames(md.14)
)

seur15.counts <- FeatureMatrix(
  fragments = frags.15,
  features = combined.peaks,
  cells = rownames(md.15)
)

seur16.counts <- FeatureMatrix(
  fragments = frags.16,
  features = combined.peaks,
  cells = rownames(md.16)
)

seur21.counts <- FeatureMatrix(
  fragments = frags.21,
  features = combined.peaks,
  cells = rownames(md.21)
)

seur22.counts <- FeatureMatrix(
  fragments = frags.22,
  features = combined.peaks,
  cells = rownames(md.22)
)

seur23.counts <- FeatureMatrix(
  fragments = frags.23,
  features = combined.peaks,
  cells = rownames(md.23)
)

seur1.counts <- FeatureMatrix(
  fragments = frags.1,
  features = combined.peaks,
  cells = rownames(md.1)
)

seur2.counts <- FeatureMatrix(
  fragments = frags.2,
  features = combined.peaks,
  cells = rownames(md.2)
)

seur3.counts <- FeatureMatrix(
  fragments = frags.3,
  features = combined.peaks,
  cells = rownames(md.3)
)

seur4.counts <- FeatureMatrix(
  fragments = frags.4,
  features = combined.peaks,
  cells = rownames(md.4)
)

seur5.counts <- FeatureMatrix(
  fragments = frags.5,
  features = combined.peaks,
  cells = rownames(md.5)
)

seur6.counts <- FeatureMatrix(
  fragments = frags.6,
  features = combined.peaks,
  cells = rownames(md.6)
)

seur7.counts <- FeatureMatrix(
  fragments = frags.7,
  features = combined.peaks,
  cells = rownames(md.7)
)

seur8.counts <- FeatureMatrix(
  fragments = frags.8,
  features = combined.peaks,
  cells = rownames(md.8)
)

seur01_assay <- CreateChromatinAssay(seur01.counts, fragments = frags.01)
seur01_ATAC <- CreateSeuratObject(seur01_assay, assay = "ATAC", meta.data = md.01, project = "CTRL1")

seur02_assay <- CreateChromatinAssay(seur02.counts, fragments = frags.02)
seur02_ATAC <- CreateSeuratObject(seur02_assay, assay = "ATAC", meta.data = md.02, project = "CTRL2")

seur03_assay <- CreateChromatinAssay(seur03.counts, fragments = frags.03)
seur03_ATAC <- CreateSeuratObject(seur03_assay, assay = "ATAC", meta.data = md.03, project = "CTRL3")

seur04_assay <- CreateChromatinAssay(seur04.counts, fragments = frags.04)
seur04_ATAC <- CreateSeuratObject(seur04_assay, assay = "ATAC", meta.data = md.04, project = "CTRL4")

seur70_assay <- CreateChromatinAssay(seur70.counts, fragments = frags.70)
seur70_ATAC <- CreateSeuratObject(seur70_assay, assay = "ATAC", meta.data = md.70, project = "C9noALSnoFTLD")

seur11_assay <- CreateChromatinAssay(seur11.counts, fragments = frags.11)
seur11_ATAC <- CreateSeuratObject(seur11_assay, assay = "ATAC", meta.data = md.11, project = "C9ALSFTLD1")

seur12_assay <- CreateChromatinAssay(seur12.counts, fragments = frags.12)
seur12_ATAC <- CreateSeuratObject(seur12_assay, assay = "ATAC", meta.data = md.12, project = "C9ALSFTLD2")

seur13_assay <- CreateChromatinAssay(seur13.counts, fragments = frags.13)
seur13_ATAC <- CreateSeuratObject(seur13_assay, assay = "ATAC", meta.data = md.13, project = "C9ALSFTLD3")

seur14_assay <- CreateChromatinAssay(seur14.counts, fragments = frags.14)
seur14_ATAC <- CreateSeuratObject(seur14_assay, assay = "ATAC", meta.data = md.14, project = "C9ALSFTLD5")

seur15_assay <- CreateChromatinAssay(seur15.counts, fragments = frags.15)
seur15_ATAC <- CreateSeuratObject(seur15_assay, assay = "ATAC", meta.data = md.15, project = "C9ALSFTLD6")

seur16_assay <- CreateChromatinAssay(seur16.counts, fragments = frags.16)
seur17_ATAC <- CreateSeuratObject(seur16_assay, assay = "ATAC", meta.data = md.16, project = "C9ALSFTLD7")

seur21_assay <- CreateChromatinAssay(seur21.counts, fragments = frags.21)
seur21_ATAC <- CreateSeuratObject(seur21_assay, assay = "ATAC", meta.data = md.21, project = "C9ALSnoFTLD1")

seur22_assay <- CreateChromatinAssay(seur22.counts, fragments = frags.22)
seur22_ATAC <- CreateSeuratObject(seur22_assay, assay = "ATAC", meta.data = md.22, project = "C9ALSnoFTLD2")

seur23_assay <- CreateChromatinAssay(seur23.counts, fragments = frags.23)
seur23_ATAC <- CreateSeuratObject(seur23_assay, assay = "ATAC", meta.data = md.23, project = "C9ALSnoFTLD3")

seur1_assay <- CreateChromatinAssay(seur1.counts, fragments = frags.1)
seur1_ATAC <- CreateSeuratObject(seur1_assay, assay = "ATAC", meta.data = md.1, project = "sALSnoFTLD1")

seur2_assay <- CreateChromatinAssay(seur2.counts, fragments = frags.2)
seur2_ATAC <- CreateSeuratObject(seur2_assay, assay = "ATAC", meta.data = md.2, project = "sALSnoFTLD2")

seur3_assay <- CreateChromatinAssay(seur3.counts, fragments = frags.3)
seur3_ATAC <- CreateSeuratObject(seur3_assay, assay = "ATAC", meta.data = md.3, project = "sALSnoFTLD3")

seur4_assay <- CreateChromatinAssay(seur4.counts, fragments = frags.4)
seur4_ATAC <- CreateSeuratObject(seur4_assay, assay = "ATAC", meta.data = md.4, project = "sALSnoFTLD4")

seur5_assay <- CreateChromatinAssay(seur5.counts, fragments = frags.5)
seur5_ATAC <- CreateSeuratObject(seur5_assay, assay = "ATAC", meta.data = md.5, project = "sALSnoFTLD5")

seur6_assay <- CreateChromatinAssay(seur6.counts, fragments = frags.6)
seur6_ATAC <- CreateSeuratObject(seur6_assay, assay = "ATAC", meta.data = md.6, project = "sALSnoFTLD6")

seur7_assay <- CreateChromatinAssay(seur7.counts, fragments = frags.7)
seur7_ATAC <- CreateSeuratObject(seur7_assay, assay = "ATAC", meta.data = md.7, project = "sALSnoFTLD7")

seur8_assay <- CreateChromatinAssay(seur8.counts, fragments = frags.8)
seur8_ATAC <- CreateSeuratObject(seur8_assay, assay = "ATAC", meta.data = md.8, project = "sALSnoFTLD8")

seur01_ATAC$sample <- 'CTRL1'
seur02_ATAC$sample <- 'CTRL2'
seur03_ATAC$sample <- 'CTRL3'
seur04_ATAC$sample <- 'CTRL4'

seur70_ATAC$sample <- 'C9noALSnoFTLD'

seur1_ATAC$sample <- 'sALSnoFTLD1'
seur2_ATAC$sample <- 'sALSnoFTLD2'
seur3_ATAC$sample <- 'sALSnoFTLD3'
seur4_ATAC$sample <- 'sALSnoFTLD4'
seur5_ATAC$sample <- 'sALSnoFTLD5'
seur6_ATAC$sample <- 'sALSnoFTLD6'
seur7_ATAC$sample <- 'sALSnoFTLD7'
seur8_ATAC$sample <- 'sALSnoFTLD8'

seur11_ATAC$sample <- 'C9ALSFTLD1'
seur12_ATAC$sample <- 'C9ALSFTLD2'
seur13_ATAC$sample <- 'C9ALSFTLD3'
seur14_ATAC$sample <- 'C9ALSFTLD5'
seur15_ATAC$sample <- 'C9ALSFTLD6'
seur16_ATAC$sample <- 'C9ALSFTLD7'

seur21_ATAC$sample <- 'C9ALSnoFTLD1'
seur22_ATAC$sample <- 'C9ALSnoFTLD2'
seur23_ATAC$sample <- 'C9ALSnoFTLD3'

seur01_ATAC$diagnosis <- 'control'
seur02_ATAC$diagnosis <- 'control'
seur03_ATAC$diagnosis <- 'control'
seur04_ATAC$diagnosis <- 'control'

seur70_ATAC$diagnosis <- 'C9noALSnoFTLD'

seur1_ATAC$diagnosis <- 'sALSnoFTLD'
seur2_ATAC$diagnosis <- 'sALSnoFTLD'
seur3_ATAC$diagnosis <- 'sALSnoFTLD'
seur4_ATAC$diagnosis <- 'sALSnoFTLD'
seur5_ATAC$diagnosis <- 'sALSnoFTLD'
seur6_ATAC$diagnosis <- 'sALSnoFTLD'
seur7_ATAC$diagnosis <- 'sALSnoFTLD'
seur8_ATAC$diagnosis <- 'sALSnoFTLD'

seur11_ATAC$diagnosis <- 'C9ALSFTLD'
seur12_ATAC$diagnosis <- 'C9ALSFTLD'
seur13_ATAC$diagnosis <- 'C9ALSFTLD'
seur14_ATAC$diagnosis <- 'C9ALSFTLD'
seur15_ATAC$diagnosis <- 'C9ALSFTLD'
seur16_ATAC$diagnosis <- 'C9ALSFTLD'

seur21_ATAC$diagnosis <- 'C9ALSnoFTLD'
seur22_ATAC$diagnosis <- 'C9ALSnoFTLD'
seur23_ATAC$diagnosis <- 'C9ALSnoFTLD'

seur01_ATAC$sex <- "female"
seur02_ATAC$sex <- "male"
seur03_ATAC$sex <- "female"
seur04_ATAC$sex <- "male"
seur70_ATAC$sex <- "male"

seur1_ATAC$sex <- "male"
seur2_ATAC$sex <- "female"
seur3_ATAC$sex <- "male"
seur4_ATAC$sex <- "female"
seur5_ATAC$sex <- "male"
seur6_ATAC$sex <- "female"
seur7_ATAC$sex <- "male"
seur8_ATAC$sex <- "male"

seur11_ATAC$sex <- "male"
seur12_ATAC$sex <- "male"
seur13_ATAC$sex <- "male"
seur14_ATAC$sex <- "female"
seur15_ATAC$sex <- "male"
seur16_ATAC$sex <- "male"
seur21_ATAC$sex <- "female"
seur22_ATAC$sex <- "female"
seur23_ATAC$sex <- "female"

system("mkdir objects")
saveRDS(seur01_ATAC,file="objects/seur01_ATAC.RDS")
saveRDS(seur02_ATAC,file="objects/seur02_ATAC.RDS")
saveRDS(seur03_ATAC,file="objects/seur03_ATAC.RDS")
saveRDS(seur04_ATAC,file="objects/seur04_ATAC.RDS")
saveRDS(seur70_ATAC,file="objects/seur70_ATAC.RDS")
saveRDS(seur1_ATAC,file="objects/seur1_ATAC.RDS")
saveRDS(seur2_ATAC,file="objects/seur2_ATAC.RDS")
saveRDS(seur3_ATAC,file="objects/seur3_ATAC.RDS")
saveRDS(seur4_ATAC,file="objects/seur4_ATAC.RDS")
saveRDS(seur5_ATAC,file="objects/seur5_ATAC.RDS")
saveRDS(seur6_ATAC,file="objects/seur6_ATAC.RDS")
saveRDS(seur7_ATAC,file="objects/seur7_ATAC.RDS")
saveRDS(seur8_ATAC,file="objects/seur8_ATAC.RDS")
saveRDS(seur11_ATAC,file="objects/seur11_ATAC.RDS")
saveRDS(seur12_ATAC,file="objects/seur12_ATAC.RDS")
saveRDS(seur13_ATAC,file="objects/seur13_ATAC.RDS")
saveRDS(seur14_ATAC,file="objects/seur14_ATAC.RDS")
saveRDS(seur15_ATAC,file="objects/seur15_ATAC.RDS")
saveRDS(seur16_ATAC,file="objects/seur16_ATAC.RDS")
saveRDS(seur21_ATAC,file="objects/seur21_ATAC.RDS")
saveRDS(seur22_ATAC,file="objects/seur22_ATAC.RDS")
saveRDS(seur23_ATAC,file="objects/seur23_ATAC.RDS")

