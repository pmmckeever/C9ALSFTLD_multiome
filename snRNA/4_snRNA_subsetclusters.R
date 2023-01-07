library(Seurat)
library(viridis)
library(future)
set.seed(1234)

plan("multisession",workers=4)
options(future.globals.maxSize = 100 * 1024^3, seed=TRUE)
table(intRNA$seurat_clusters)

intRNA<-readRDS(file="objects/intRNA_SCT.RDS")

pdf("figures/DimPlot_celltypes.pdf", height=7, width=10)
DimPlot(
  object = intRNA,
  group.by="seurat_clusters",label=T,
  reduction="umap")
dev.off()

#check Violin plots for canonical cell type markers
pdf("figures/SYT1_snRNA.pdf", height=5, width=12)
VlnPlot(object = intRNA,features = c("SYT1") , pt.size=0) + NoLegend()
dev.off()

pdf("figures/NRGN_snRNA.pdf", height=5, width=12)
VlnPlot(object = intRNA,features = c("NRGN") , pt.size=0) + NoLegend() 
dev.off()

pdf("figures/SLC17A7_snRNA.pdf", height=5, width=12)
VlnPlot(object = intRNA,features = c("SLC17A7") , pt.size=0) + NoLegend() 
dev.off()

pdf("figures/SATB2_snRNA.pdf", height=5, width=12)
VlnPlot(object = intRNA,features = c("SLC17A7") , pt.size=0) + NoLegend() 
dev.off()

pdf("figures/GAD1_snRNA.pdf", height=5, width=12)
VlnPlot(object = intRNA,features = c("GAD1") , pt.size=0) + NoLegend()
dev.off()

pdf("figures/GAD2_snRNA.pdf", height=5, width=12)
VlnPlot(object = intRNA,features = c("GAD2") , pt.size=0) + NoLegend()
dev.off()

pdf("figures/AQP4_snRNA.pdf", height=5, width=12)
VlnPlot(object = intRNA,features = c("AQP4") , pt.size=0) + NoLegend() 
dev.off()

pdf("figures/SLC1A2_snRNA.pdf", height=5, width=12)
VlnPlot(object = intRNA,features = c("SLC1A2") , pt.size=0) + NoLegend() 
dev.off()

pdf("figures/OPALIN_snRNA.pdf", height=5, width=12)
VlnPlot(object = intRNA,features = c("OPALIN") , pt.size=0) + NoLegend() 
dev.off()

pdf("figures/PLP1_snRNA.pdf", height=5, width=12)
VlnPlot(object = intRNA,features = c("PLP1") , pt.size=0) + NoLegend()
dev.off()

pdf("figures/PDGFRA_snRNA.pdf", height=5, width=12)
VlnPlot(object = intRNA,features = c("PDGFRA") , pt.size=0) + NoLegend()
dev.off()

pdf("figures/VCAN_snRNA.pdf", height=5, width=12)
VlnPlot(object = intRNA,features = c("VCAN") , pt.size=0) + NoLegend() 
dev.off()

pdf("figures/CD74_snRNA.pdf", height=5, width=12)
VlnPlot(object = intRNA,features = c("CD74") , pt.size=0) + NoLegend() 
dev.off()

pdf("figures/CSF1R_snRNA.pdf", height=5, width=12)
VlnPlot(object = intRNA,features = c("CSF1R") , pt.size=0) + NoLegend() 
dev.off()

pdf("figures/FLT1_snRNA.pdf", height=5, width=12)
VlnPlot(object = intRNA,features = c("FLT1") , pt.size=0) + NoLegend()
dev.off()

pdf("figures/CLDN5_snRNA.pdf", height=5, width=12)
VlnPlot(object = intRNA,features = c("CLDN5") , pt.size=0) + NoLegend() 
dev.off()

#remove suspected double cluster
intRNA<-subset(intRNA, idents=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","15","16","17",
                                "18","19","20","21","22","23","24"))

DefaultAssay(object = intRNA) <- "integrated"
#intRNA <- ScaleData(object = intRNA, verbose = FALSE) 
#intRNA <- FindVariableFeatures(intRNA)
intRNA <- RunPCA(object = intRNA, npcs = 30, verbose = FALSE)
intRNA <- RunTSNE(object = intRNA, dims=1:30, perplexity=30)
intRNA <- RunUMAP(object = intRNA, reduction = "pca", dims = 1:30)
intRNA <- FindNeighbors(intRNA, reduction="pca", dims=1:30)
intRNA <- FindClusters(intRNA, resolution=0.8)
Idents(intRNA)<-'seurat_clusters'
table(intRNA$seurat_clusters)

saveRDS(intRNA,file="objects/intRNA_noannotation.RDS")