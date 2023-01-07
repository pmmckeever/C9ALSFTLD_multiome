library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(GenomeInfoDb)
library(patchwork)
library(ggplot2)
set.seed(1234)

intATAC <-readRDS(file="objects/ATAC_nomults2.RDS")

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- 'hg38'
DefaultAssay(intATAC) <- "ATAC"
Annotation(intATAC) <- annotations

# trim chromosomes
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- as.logical(seqnames(granges(intATAC)) %in% main.chroms)
intATAC <- intATAC[keep.peaks, ]

gene.activities <- GeneActivity(intATAC)
intATAC[['RNA']] <- CreateAssayObject(counts = gene.activities)

DefaultAssay(intATAC)<-"RNA"
intATAC <- NormalizeData(
  object = intATAC,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median
(intATAC$nCount_RNA)
)
saveRDS(intATAC,file="objects/ATAC_RNAact.RDS")

intRNA<-readRDS(file="objects/Velmeshev_CTRLonly_PFC_intRNA.RDS")
DefaultAssay(intATAC)<-"RNA"
DefaultAssay(intRNA)<-"integrated"

intRNA<-FindVariableFeatures(intRNA, nfeatures=3000)
intATAC<-FindVariableFeatures(intATAC, nfeatures=3000)

transfer.anchors <- FindTransferAnchors(
  reference = intRNA,
  query = intATAC,
  dims = 1:30,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = intRNA$celltype,
  weight.reduction = intATAC[['harmony']],
  dims = 2:30
)

intATAC <- AddMetaData(object = intATAC, metadata = predicted.labels)

Idents(intATAC)<-"predicted.id"

intATAC$cellsubtype<-Idents(intATAC)

table(intATAC$cellsubtype)

intATAC2 <- intATAC[,intATAC$prediction.score.max >= 0.5]
table(intATAC2$cellsubtype)

saveRDS(intATAC2,file="objects/ATAC_labelxfer.RDS")

Idents(intATAC)<-"predicted.labels"
pdf("figures/DimPlot_ATAC_samples_labelxfer.pdf", height=7, width=11)
DimPlot(
  object = intATAC,
  group.by="sample",
  reduction="umap")
dev.off()

pdf("figures/DimPlot_ATAC_group_labelxfer.pdf", height=7, width=9)
DimPlot(
  object = intATAC,
  group.by="diagnosis",
  reduction="umap")
dev.off()

pdf("figures/DimPlot_ATAC_celltypes_labelxfer.pdf", height=7, width=9)
DimPlot(
  object = intATAC,
  group.by="predicted.id",
  reduction="umap",label=T)
dev.off()