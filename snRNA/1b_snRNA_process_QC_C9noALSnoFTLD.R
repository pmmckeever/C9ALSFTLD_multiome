library(Seurat)
library(org.Hs.eg.db)
library(DoubletFinder)
set.seed(1234)

#CTRL matrix location
input_from_10x70 <- "/ALS_snRNA/C9noALSnoFTLD_snRNA/outs/filtered_feature_bc_matrix"

#C9noALSnoFTLD cases read-in
seur70_RNA <- CreateSeuratObject(counts=Read10X(input_from_10x70), project="C9noALSnoFTLD",
                          min.cells=1,min.features=1)
show(seur70_RNA)
seur70_RNA[["percent.mt"]] <- PercentageFeatureSet(seur70_RNA, pattern = "^MT-")

mito_gene_identifier <- "^MT-"
mads_thresh <- 3
hard_thresh <- 50
seur70_RNA <- PercentageFeatureSet(seur70_RNA, pattern = "^MT-", col.name = "percent.mt")
mito_thresh <- median(seur70_RNA$percent.mt) + mad(seur70_RNA$percent.mt) * mads_thresh
drop_mito <- seur70_RNA$percent.mt > mito_thresh | seur70_RNA$percent.mt > hard_thresh
seur70_RNA_mitoBI <- seur70_RNA[,!drop_mito]
show(seur70_RNA_mitoBI)

seur70_RNA_mitoBI <- subset(seur70_RNA_mitoBI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)

seur70_RNA_mitoBI$sample <- 'C9noALSnoFTLD'

seur70_RNA_mitoBI<-NormalizeData(seur70_RNA_mitoBI)

seur70_RNA_mitoBI<-FindVariableFeatures(seur70_RNA_mitoBI)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

seur70_RNA_mitoBI <- CellCycleScoring(seur70_RNA_mitoBI, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DefaultAssay(object = seur70_RNA_mitoBI) <- "RNA"
seur70_RNA_mitoBI <- ScaleData(object = seur70_RNA_mitoBI, verbose = FALSE)
seur70_RNA_mitoBI <- FindVariableFeatures(seur70_RNA_mitoBI)
seur70_RNA_mitoBI <- RunPCA(object = seur70_RNA_mitoBI, npcs = 50, verbose = FALSE)
seur70_RNA_mitoBI <- RunUMAP(object = seur70_RNA_mitoBI, reduction = "pca", dims = 1:50)
seur70_RNA_mitoBI <- FindNeighbors(seur70_RNA_mitoBI,reduction="pca",dims=1:50)
seur70_RNA_mitoBI <- FindClusters(seur70_RNA_mitoBI,resolution=0.8)

## pK Identification (no ground-truth)
sweep.res.list_seur70_RNA_mitoBI <- paramSweep_v3(seur70_RNA_mitoBI, PCs = 1:50, sct = FALSE)
sweep.stats_seur70_RNA_mitoBI <- summarizeSweep(sweep.res.list_seur70_RNA_mitoBI, GT = FALSE)
bcmvn_seur70_RNA_mitoBI <- find.pK(sweep.stats_seur70_RNA_mitoBI)
bcmvn_seur70_RNA_mitoBI

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seur70_RNA_mitoBI@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.05*nrow(seur70_RNA_mitoBI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seur70_RNA_mitoBI <- doubletFinder_v3(seur70_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.18, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
#seur70_RNA_mitoBI <- doubletFinder_v3(seur70_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
DimPlot(seur70_RNA_mitoBI, pt.size=0.1, group.by="DF.classifications_0.25_0.18_234")

## paste column to new place and delete old one
#seur70_RNA_mitoBI$sample_doublets <- seur70_RNA_mitoBI$DF.classifications_0.25_0.18_234
#seur70_RNA_mitoBI$DF.classifications_0.25_0.18_234 <- NULL

#remove doublets and shrink object
seur70_RNA_singlet<-subset(x=seur70_RNA_mitoBI, subset = DF.classifications_0.25_0.18_234 == "Singlet")
seur70_RNA_singlet_diet<-DietSeurat(seur70_RNA_singlet)
seur70_RNA_singlet_diet$DF.classifications_0.25_0.18_234<-NULL
seur70_RNA_singlet_diet$pANN_0.25_0.18_234<-NULL
seur70_RNA_singlet_diet$seurat_clusters<-NULL
seur70_RNA_singlet_diet$RNA_snn_res.0.8<-NULL
saveRDS(seur70_RNA_singlet_diet,file="objects/seur70_RNA.RDS")