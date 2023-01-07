library(Signac)
library(Seurat)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(dplyr)
library(tidyr)
library(ggplot2)
set.seed(1234)

intATAC <- readRDS(file="objects/ATAC_labelxfer.RDS")

DefaultAssay(intATAC) <- "ATAC"
Idents(intATAC) <- "cellsubtype"

peaks <-CallPeaks(intATAC,macs2.path="/ENV/lib/python3.9/site-packages/MACS2",group.by="cellsubtype")
peaks2 <- keepStandardChromosomes(peaks, pruning.mode="coarse")
peaks2 <- subsetByOverlaps(x=peaks2, ranges=blacklist_hg38_unified, invert=T)
macs2_counts <- FeatureMatrix(fragments=Fragments(intATAC), features=peaks2, cells=colnames(intATAC))

#save in case crash
saveRDS(macs2_counts, file="objects/macs2counts.RDS")
saveRDS(peaks2, file="objects/ATAC_MACSpeaks.RDS")

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

intATAC[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = Fragments(intATAC),
  annotation = annotations
)

saveRDS(intATAC,file="objects/ATAC_MACS2.RDS")

pfm <- getMatrixSet(x = JASPAR2020,opts = list(species = 9606, all_versions = FALSE))
intATAC <- AddMotifs(intATAC,genome=BSgenome.Hsapiens.UCSC.hg38, pfm=pfm)
gc()
intATAC <- RunChromVAR(object = intATAC,genome = BSgenome.Hsapiens.UCSC.hg38)
intATAC

intATAC<-RegionStats(intATAC,genome=BSgenome.Hsapiens.UCSC.hg38)
intATAC

#processing step. not used for visualization... replace with Cicero linked peaks
intATAC<-LinkPeaks(intATAC,peak.assay="peaks",expression.assay="RNA")
intATAC

saveRDS(intATAC,file="objects/ATAC_SignacPeaks.RDS")
