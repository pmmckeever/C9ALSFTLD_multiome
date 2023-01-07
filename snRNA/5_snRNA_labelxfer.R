library(Seurat)
set.seed(1234)

## read in Velmeshev et al data and save control data into new object
other.data <- Read10X(data.dir = "/other_data/Velmeshev_et_al/rawMatrix")
meta.df <- read.delim(file.path("/other_data/Velmeshev_et_al/rawMatrix", "meta.txt"))
meta.df2<-meta.df[,-1]
rownames(meta.df2)<-meta.df[,1]
velm <- CreateSeuratObject(other.data, meta.data = meta.df2)
Idents(velm)<-"diagnosis"
control<-subset(velm,idents="Control")
control_PFC<-subset(velm,idents="PFC")
Idents(velm)<-"region"
control_PFC<-subset(velm,idents="PFC")
saveRDS(control_PFC,file="/other_data/Velmeshev_et_al/Velmeshev_control_PFC_RNAonly.RDS")

DefaultAssay(controlPFC)<-"RNA"

#split file
PFC.list<-SplitObject(controlPFC, split.by="sample")

# normalize and identify variable features for each dataset independently
PFC.list <- lapply(X = PFC.list, FUN = SCTransform)

#integrate based on features (SCT)
features <- SelectIntegrationFeatures(object.list = PFC.list, nfeatures=3000)
PFC.list <- PrepSCTIntegration(object.list = PFC.list, anchor.features=features)

#
anchors <- FindIntegrationAnchors(object.list = PFC.list, normalization.method="SCT", reference = c(5, 12), anchor.features=features, k.anchor=20)
intPFC <- IntegrateData(anchorset = anchors, normalization.method="SCT")
DefaultAssay(intPFC) <- "integrated"

#dimensionality reduction and clustering
intPFC <- RunPCA(object = intPFC, npcs = 50, verbose = FALSE)
intPFC <- RunTSNE(object = intPFC, dims=1:50, perplexity=30)
intPFC <- RunUMAP(object = intPFC, reduction = "pca", dims = 1:50, min.dist = 0.1)
intPFC <- FindNeighbors(intPFC, reduction="pca", dims=1:50)
intPFC <- FindClusters(intPFC, resolution=0.6)
saveRDS(intPFC,file="objects/Velmeshev_CTRLonly_PFC_intRNA.RDS")

intRNA <- readRDS(file="objects/intRNA_noannotation.RDS")
DefaultAssay(intRNA)<-"integrated"
DefaultAssay(intPFC)<-"integrated"

transfer.anchors <- FindTransferAnchors(
  reference = intPFC,
  query = intRNA,
  dims = 1:30,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = intPFC$cluster,
  weight.reduction = 'cca',
  dims = 1:30
)

intRNA <- AddMetaData(object = intRNA, metadata = predicted.labels)

Idents(intRNA)<-"predicted.id"

intRNA$cellsubtype<-Idents(intRNA)

table(intRNA$cellsubtype)

saveRDS(intRNA,file="objects/intRNA_labelxfer_nocutoff.RDS")

intRNA2 <- intRNA[,intRNA$prediction.score.max >= 0.5]
table(intRNA$cellsubtype)

saveRDS(intRNA2,file="objects/intRNA_labelxfer.RDS")

Idents(intRNA)<-"predicted.labels"
pdf("figures/DimPlot_RNA_cellsubtypes_labelxfer.pdf", height=7, width=11)
DimPlot(
  object = intRNA,
  group.by="cellsubtype",
  reduction="umap")
dev.off()