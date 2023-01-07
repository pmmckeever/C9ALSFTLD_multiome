library(Seurat)
set.seed(1234)

merger <- readRDS(file="objects/unint_RNA.RDS")
DefaultAssay(merger) <- "RNA"
merger$diagnosis<-merger$group

#split file
merger.list <- SplitObject(merger, split.by="sample")

# normalize and identify variable features for each dataset independently
merger.list <- lapply(X = merger.list, FUN = SCTransform), vars.to.regress = c("percent.mt"))

#integrate based on features (SCT)
features <- SelectIntegrationFeatures(object.list = merger.list, nfeatures=3000)
merger.list <- PrepSCTIntegration(object.list = merger.list, anchor.features=features)

#change k.anchor to increase overlap
anchors <- FindIntegrationAnchors(object.list = merger.list, normalization.method="SCT", reference = c(1, 2), anchor.features=features, k.anchor=20)
merger.integrated <- IntegrateData(anchorset = anchors, normalization.method="SCT")
DefaultAssay(merger.integrated) <- "integrated"

#dimensionality reduction and clustering
merger.integrated <- RunPCA(object = merger.integrated, npcs = 50, verbose = FALSE)
merger.integrated <- RunTSNE(object = merger.integrated, dims=1:50, perplexity=30)
merger.integrated <- RunUMAP(object = merger.integrated, reduction = "pca", dims = 1:50, min.dist = 0.1)
merger.integrated <- FindNeighbors(merger.integrated, reduction="pca",dims=1:50)
merger.integrated <- FindClusters(merger.integrated, resolution=0.6)
DefaultAssay(merger.integrated)<-"RNA"
saveRDS(merger.integrated,file="objects/intRNA_SCT.RDS")

#visualize
DefaultAssay(merger.integrated) <- "SCT"

pdf("figures/DimPlot_merger.integrated_batch.pdf", height=7, width=11)
DimPlot(
  object = merger.integrated,
  group.by="sample",
  reduction="umap")
dev.off()

pdf("figures/DimPlot_merger.integrated_group.pdf", height=7, width=9)
DimPlot(
  object = merger.integrated,
  group.by="diagnosis",
  reduction="umap")
dev.off()

pdf("figures/DimPlot_merger.integrated_sex.pdf", height=7, width=9)
DimPlot(
  object = merger.integrated,
  group.by="sex",
  reduction="umap")
dev.off()

pdf("figures/DimPlot_merger.integrated_celltypes.pdf", height=7, width=10)
DimPlot(
  object = merger.integrated,
  group.by="seurat_clusters",
  reduction="umap")
dev.off()
