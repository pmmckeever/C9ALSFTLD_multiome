library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(patchwork)
set.seed(1234)

seur01_ATAC<-readRDS(file="objects/seur01_ATAC.RDS")
seur02_ATAC<-readRDS(file="objects/seur02_ATAC.RDS")
seur03_ATAC<-readRDS(file="objects/seur03_ATAC.RDS")
seur04_ATAC<-readRDS(file="objects/seur04_ATAC.RDS")
seur70_ATAC<-readRDS(file="objects/seur70_ATAC.RDS")
seur1_ATAC<-readRDS(file="objects/seur1_ATAC.RDS")
seur2_ATAC<-readRDS(file="objects/seur2_ATAC.RDS")
seur3_ATAC<-readRDS(file="objects/seur3_ATAC.RDS")
seur4_ATAC<-readRDS(file="objects/seur4_ATAC.RDS")
seur5_ATAC<-readRDS(file="objects/seur5_ATAC.RDS")
seur6_ATAC<-readRDS(file="objects/seur6_ATAC.RDS")
seur7_ATAC<-readRDS(file="objects/seur7_ATAC.RDS")
seur8_ATAC<-readRDS(file="objects/seur8_ATAC.RDS")
seur11_ATAC<-readRDS(file="objects/seur11_ATAC.RDS")
seur12_ATAC<-readRDS(file="objects/seur12_ATAC.RDS")
seur13_ATAC<-readRDS(file="objects/seur13_ATAC.RDS")
seur15_ATAC<-readRDS(file="objects/seur15_ATAC.RDS")
seur16_ATAC<-readRDS(file="objects/seur16_ATAC.RDS")
seur17_ATAC<-readRDS(file="objects/seur17_ATAC.RDS")
seur21_ATAC<-readRDS(file="objects/seur21_ATAC.RDS")
seur22_ATAC<-readRDS(file="objects/seur22_ATAC.RDS")
seur23_ATAC<-readRDS(file="objects/seur23_ATAC.RDS")

#create object without integration first
unint_ATAC <- merge(
  x = seur01_ATAC,
  y = list(seur02_ATAC, seur03_ATAC, seur04_ATAC, seur70_ATAC,
    seur1_ATAC, seur2_ATAC, seur3_ATAC, seur4_ATAC, seur5_ATAC, seur6_ATAC, seur7_ATAC, seur8_ATAC,
    seur11_ATAC, seur12_ATAC, seur13_ATAC, seur15_ATAC, seur16_ATAC, seur17_ATAC,
    seur21_ATAC, seur22_ATAC, seur23_ATAC), add.cell.ids=c("CTRL1","CTRL2","CTRL3","CTRL4","C9noALSnoFTLD",
                    "sALSnoFTLD1","sALSnoFTLD2","sALSnoFTLD3","sALSnoFTLD4","sALSnoFTLD5","sALSnoFTLD6","sALSnoFTLD7","sALSnoFTLD8",
                    "C9ALSFTLD1","C9ALSFTLD2","C9ALSFTLD3","C9ALSFTLD5","C9ALSFTLD6","C9ALSFTLD7",
                    "C9ALSnoFTLD1","C9ALSnoFTLD2","C9ALSnoFTLD3"), project="snATAC")

rm(seurC1_ATAC,seurC2_ATAC,seurC3_ATAC,seurC4_ATAC,seur70_ATAC,seur1_ATAC,seur2_ATAC,seur3_ATAC,seur4_ATAC,
  seur5_ATAC,seur6_ATAC,seur7_ATAC,seur8_ATAC,seur11_ATAC,seur12_ATAC,seur13_ATAC,seur15_ATAC,seur16_ATAC,
  seur17_ATAC,seur21_ATAC,seur22_ATAC,seur23_ATAC)
gc()

unint_ATAC[["ATAC"]]
unint_ATAC <- FindVariableFeatures(unint_ATAC)
unint_ATAC <- RunTFIDF(unint_ATAC)
unint_ATAC <- FindTopFeatures(unint_ATAC, min.cutoff = 50)
unint_ATAC <- RunSVD(unint_ATAC)
unint_ATAC <- RunUMAP(unint_ATAC, reduction = 'lsi', dims = 2:30)

#inspect plots for batch effects
system("mkdir figures")
pdf("figures/DimPlot_noQC_noint_ATAC.pdf")
DimPlot(unint_ATAC,pt.size=0.1,group.by="sample")
dev.off()
pdf("figures/DimPlot_noQC_noint_ATAC.pdf")
DimPlot(unint_ATAC,pt.size=0.1,group.by="diagnosis")
dev.off()
gc()

DefaultAssay(unint_ATAC) <- "ATAC"

granges(unint_ATAC)
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(unint_ATAC) <- annotations

# compute nucleosome signal score per cell
unint_ATAC <- NucleosomeSignal(object = unint_ATAC)

# compute TSS enrichment score per cell
unint_ATAC <- TSSEnrichment(object = unint_ATAC, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
unint_ATAC$pct_reads_in_peaks <- unint_ATAC$peak_region_fragments / unint_ATAC$passed_filters * 100
unint_ATAC$blacklist_ratio <- unint_ATAC$blacklist_region_fragments / unint_ATAC$peak_region_fragments
unint_ATAC$high.tss <- ifelse(unint_ATAC$TSS.enrichment > 2, 'High', 'Low')
unint_ATAC$nucleosome_group <- ifelse(unint_ATAC$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

saveRDS(unint_ATAC,file="objects/unint_ATAC.RDS")

#more QC plots
pdf("figures/TSSplot.pdf")
TSSPlot(unint_ATAC, group.by = 'high.tss') + NoLegend()
dev.off()

pdf("figures/FragHisto.pdf")
FragmentHistogram(object = unint_ATAC, group.by = 'nucleosome_group')
dev.off()

pdf("figures/VlnPlot_QC1.pdf", height=7, width=14)
VlnPlot(
  object = unint_ATAC,
  features = c('pct_reads_in_peaks'),
  pt.size = 0.1
)
dev.off()

pdf("figures/VlnPlot_QC2.pdf", height=7, width=14)
VlnPlot(
  object = unint_ATAC,
  features = c('peak_region_fragments'),
  pt.size = 0.1
)
dev.off()

pdf("figures/VlnPlot_QC3_3.pdf", height=7, width=14)
VlnPlot(
  object = unint_ATAC,
  features = c('TSS.enrichment'),
  pt.size = 0.1
)
dev.off()

pdf("figures/VlnPlot_QC4_3.pdf", height=7, width=14)
VlnPlot(
  object = unint_ATAC,
  features = c('nucleosome_signal'),
  pt.size = 0.1
)
dev.off()

