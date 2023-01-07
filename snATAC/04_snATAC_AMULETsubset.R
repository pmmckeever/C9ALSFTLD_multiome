library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(GenomeInfoDb)
library(patchwork)
library(ggplot2)
library(harmony)
library(future)
set.seed(1234)

plan("multisession",workers=32)
options(future.globals.maxSize = 185 * 1024^3, seed=TRUE)


intATAC <- readRDS("objects/hmint_ATAC.RDS")

multiplets <- read.table("All_MultipletCellIds.txt") %>% t() %>% as.vector()
head(multiplets)
cells <- Cells(intATAC)
multiplet.data <- rep("singlet",length(cells))
names(multiplet.data) <- cells
multiplet.data[multiplets] <- "multiplet"

intATAC <- AddMetaData(intATAC, multiplet.data, col.name = "multiplet")

#inspect AMULET plot
pdf("figures/DimPlot_ATAC_multiplets.pdf", height=7, width=9)
DimPlot(intATAC, group.by = "multiplet", pt.size = 1, 
        cols = c("singlet" = "grey", "multiplet" = "purple"), 
        order = c("multiplet", "singlet"))
dev.off()

intATAC2 <- subset(intATAC, multiplet == "singlet")
intATAC2$multiplet <- NULL
rm(intATAC)

granges(intATAC)
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

# add the gene information to the object
Annotation(intATAC) <- annotations

# trim chromosomes
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- as.logical(seqnames(granges(intATAC)) %in% main.chroms)
intATAC2 <- intATAC[keep.peaks, ]
rm(intATAC)
gc()
intATAC2 <- DietSeurat(intATAC2)
intATAC2[["ATAC"]]
intATAC2 <- FindVariableFeatures(intATAC2)
intATAC2 <- RunTFIDF(intATAC2)
intATAC2 <- FindTopFeatures(intATAC2, min.cutoff = 50)
intATAC2 <- RunSVD(intATAC2)
intATAC2 <- RunUMAP(intATAC2, reduction = 'lsi', dims = 2:30)
DimPlot(intATAC2,group.by="sample")
intATAC2
gc()


intATAC2 <- RunHarmony(
object = intATAC2,
group.by.vars = 'sample',
reduction = 'lsi',
assay.use = 'ATAC',
project.dim = FALSE
)

intATAC2 <- RunUMAP(intATAC2, reduction = "harmony", dims = 2:30)
intATAC2 <- FindNeighbors(object = intATAC2, reduction = 'harmony', dims = 2:30)
intATAC2 <- FindClusters(object = intATAC2, verbose = FALSE, algorithm = 3)

DimPlot(intATAC2,group.by="sample",reduction="umap")

pdf("figures/DimPlot_ATAC_samples_recluster.pdf", height=7, width=9)
DimPlot(
  object = intATAC2,
  group.by="sample",
  reduction="umap")
dev.off()

pdf("figures/DimPlot_ATAC_diagnosis_recluster.pdf", height=7, width=9)
DimPlot(
  object = intATAC2,
  group.by="diagnosis",
  reduction="umap")
dev.off()

pdf("figures/DimPlot_ATAC_reclusters.pdf", height=7, width=9)
DimPlot(
  object = intATAC2,
  group.by="seurat_clusters",
  reduction="umap",label=T)
dev.off()

saveRDS(intATAC2,file="objects/ATAC_nomults.RDS")

#remove 26 and 27, suspected doublets not detected by AMULET
intATAC<-subset(intATAC2, idents=c("0","1","2","3","4","5","6","7","8","9","10",
                                    "11","12","13","14","15","16","17","18","19","20",
                                    "21","22","23","24","25","28","29"))


intATAC <- DietSeurat(intATAC)
intATAC[["ATAC"]]
intATAC <- FindVariableFeatures(intATAC)
intATAC <- RunTFIDF(intATAC)
intATAC <- FindTopFeatures(intATAC, min.cutoff = 50)
intATAC <- RunSVD(intATAC)
intATAC <- RunUMAP(intATAC, reduction = 'lsi', dims = 2:30)
DimPlot(intATAC,group.by="sample")
intATAC
gc()

intATAC <- RunHarmony(
object = intATAC,
group.by.vars = 'sample',
reduction = 'lsi',
assay.use = 'ATAC',
project.dim = FALSE
)

intATAC <- RunUMAP(intATAC, reduction = "harmony", dims = 2:30)
intATAC <- FindNeighbors(object = intATAC, reduction = 'harmony', dims = 2:30)
intATAC <- FindClusters(object = intATAC, verbose = FALSE, algorithm = 3)

table(intATAC$seurat_clusters)

pdf("figures/DimPlot_ATAC_samples_recluster2.pdf", height=7, width=9)
DimPlot(
  object = intATAC,
  group.by="sample",
  reduction="umap")
dev.off()

pdf("figures/DimPlot_ATAC_diagnosis_recluster2.pdf", height=7, width=9)
DimPlot(
  object = intATAC,
  group.by="diagnosis",
  reduction="umap")
dev.off()

pdf("figures/DimPlot_ATAC_reclusters2.pdf", height=7, width=9)
DimPlot(
  object = intATAC,
  group.by="seurat_clusters",
  reduction="umap",label=T)
dev.off()

saveRDS(intATAC,file="objects/ATAC_nomults2.RDS")