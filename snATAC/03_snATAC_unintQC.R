library(Seurat)
library(Signac)
library(harmony)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(GenomeInfoDb)
library(patchwork)
set.seed(1234)

unint<-readRDS(file="objects/unint_ATAC.RDS")

unint_QC <- subset(
  x = unint,
  subset = peak_region_fragments > 2000 &
    peak_region_fragments < 25000 &
    pct_reads_in_peaks > 15 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    TSS.enrichment < 20
)
unint_QC

rm(unint)
gc()

unint_QC[["ATAC"]]
unint_QC <- FindVariableFeatures(unint_QC)
unint_QC <- RunTFIDF(unint_QC)
unint_QC <- FindTopFeatures(unint_QC, min.cutoff = 50)
unint_QC <- RunSVD(unint_QC)
unint_QC <- RunUMAP(unint_QC, reduction = 'lsi', dims = 2:30)

unint_QC <- RunHarmony(
  object = unint_QC,
  group.by.vars = 'sample',
  reduction = 'lsi',
  assay.use = 'ATAC',
  project.dim = FALSE
)

unint_QC <- RunUMAP(unint_QC, reduction = "harmony", dims = 2:30)
unint_QC <- FindNeighbors(object = unint_QC, reduction = 'harmony', dims = 2:30)
unint_QC <- FindClusters(object = unint_QC, verbose = FALSE, algorithm = 3)

saveRDS(unint_QC,file="objects/hmint_ATAC.RDS")

#inspect Harmony plots

pdf("figures/DimPlot_samples.pdf", height=7, width=9)
DimPlot(
  object = unint_QC,
  group.by="sample")
dev.off()

pdf("figures/DimPlot_group.pdf", height=7, width=9)
DimPlot(
  object = unint_QC,
  group.by="diagnosis")
dev.off()

pdf("figures/DimPlot_clusters.pdf", height=7, width=9)
DimPlot(
  object = unint_QC,
  group.by="seurat_clusters",
  reduction="umap",label=T)
dev.off()
