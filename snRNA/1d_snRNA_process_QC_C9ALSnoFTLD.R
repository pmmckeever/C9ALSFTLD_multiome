library(Seurat)
library(ggplot2)
library(viridis)
library(patchwork)
library(scClustViz)
library(org.Hs.eg.db)
library(DoubletFinder)
library(future)
set.seed(1234)
plan("multicore",workers=48)
options(future.globals.maxSize = 187 * 1024^3, seed=TRUE)

#C9ALSnoFTLD matrix location - HOME SERVER
input_from_10x21 <- "/ALS_snRNA/C9ALSnoFTLD1_snRNA/outs/filtered_feature_bc_matrix"
input_from_10x22 <- "/ALS_snRNA/C9ALSnoFTLD2_snRNA/outs/filtered_feature_bc_matrix"
input_from_10x23 <- "/ALS_snRNA/C9ALSnoFTLD3_snRNA/outs/filtered_feature_bc_matrix"

#C9ALSnoFTLD cases read-in
seur21_RNA <- CreateSeuratObject(counts=Read10X(input_from_10x21), project="C9ALSnoFTLD1",
                          min.cells=1,min.features=1)
seur22_RNA <- CreateSeuratObject(counts=Read10X(input_from_10x22), project="C9ALSnoFTLD2",
                          min.cells=1,min.features=1)
seur23_RNA <- CreateSeuratObject(counts=Read10X(input_from_10x23), project="C9ALSnoFTLD3",
                            min.cells=1,min.features=1)

show(seur21_RNA)
show(seur22_RNA)
show(seur23_RNA)
show(seur24_RNA)

seur21_RNA[["percent.mt"]] <- PercentageFeatureSet(seur21_RNA, pattern = "^MT-")
seur22_RNA[["percent.mt"]] <- PercentageFeatureSet(seur22_RNA, pattern = "^MT-")
seur23_RNA[["percent.mt"]] <- PercentageFeatureSet(seur23_RNA, pattern = "^MT-")

seur21_RNA[["percent.rps"]] <- PercentageFeatureSet(seur21_RNA, pattern = "^RPS")
seur22_RNA[["percent.rps"]] <- PercentageFeatureSet(seur22_RNA, pattern = "^RPS")
seur23_RNA[["percent.rps"]] <- PercentageFeatureSet(seur23_RNA, pattern = "^RPS")

seur21_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur21_RNA, pattern = "^RPL")
seur22_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur22_RNA, pattern = "^RPL")
seur23_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur23_RNA, pattern = "^RPL")

mito_gene_identifier <- "^MT-"
mads_thresh <- 3
hard_thresh <- 50
seur21_RNA <- PercentageFeatureSet(seur21_RNA, pattern = "^MT-", col.name = "percent.mt")
mito_thresh <- median(seur21_RNA$percent.mt) + mad(seur21_RNA$percent.mt) * mads_thresh
drop_mito <- seur21_RNA$percent.mt > mito_thresh | seur21_RNA$percent.mt > hard_thresh
seur21_RNA_mitoBI <- seur21_RNA[,!drop_mito]
show(seur21_RNA_mitoBI)
seur22_RNA <- PercentageFeatureSet(seur22_RNA, pattern = "^MT-", col.name = "percent.mt")
mito_thresh <- median(seur22_RNA$percent.mt) + mad(seur22_RNA$percent.mt) * mads_thresh
drop_mito <- seur22_RNA$percent.mt > mito_thresh | seur22_RNA$percent.mt > hard_thresh
seur22_RNA_mitoBI <- seur22_RNA[,!drop_mito]
show(seur22_RNA_mitoBI)
seur23_RNA <- PercentageFeatureSet(seur23_RNA, pattern = "^MT-", col.name = "percent.mt")
mito_thresh <- median(seur23_RNA$percent.mt) + mad(seur23_RNA$percent.mt) * mads_thresh
drop_mito <- seur23_RNA$percent.mt > mito_thresh | seur23_RNA$percent.mt > hard_thresh
seur23_RNA_mitoBI <- seur23_RNA[,!drop_mito]
show(seur23_RNA_mitoBI)

seur21_RNA_mitoBI <- subset(seur21_RNA_mitoBI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)
seur22_RNA_mitoBI <- subset(seur22_RNA_mitoBI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)
seur23_RNA_mitoBI <- subset(seur23_RNA_mitoBI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)

show(seur21_RNA_mitoBI)
show(seur22_RNA_mitoBI)    
show(seur23_RNA_mitoBI)

seur21_RNA_mitoBI$sample <- 'C9ALSnoFTLD1'
seur22_RNA_mitoBI$sample <- 'C9ALSnoFTLD2'
seur23_RNA_mitoBI$sample <- 'C9ALSnoFTLD3'

seur21_RNA_mitoBI<-NormalizeData(seur21_RNA_mitoBI)
seur22_RNA_mitoBI<-NormalizeData(seur22_RNA_mitoBI)    
seur23_RNA_mitoBI<-NormalizeData(seur23_RNA_mitoBI)

seur21_RNA_mitoBI<-FindVariableFeatures(seur21_RNA_mitoBI)
seur22_RNA_mitoBI<-FindVariableFeatures(seur22_RNA_mitoBI)    
seur23_RNA_mitoBI<-FindVariableFeatures(seur23_RNA_mitoBI)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

seur21_RNA_mitoBI <- CellCycleScoring(seur21_RNA_mitoBI, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seur22_RNA_mitoBI <- CellCycleScoring(seur22_RNA_mitoBI, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seur23_RNA_mitoBI <- CellCycleScoring(seur23_RNA_mitoBI, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DefaultAssay(object = seur21_RNA_mitoBI) <- "RNA"
seur21_RNA_mitoBI <- ScaleData(object = seur21_RNA_mitoBI, verbose = FALSE)
seur21_RNA_mitoBI <- FindVariableFeatures(seur21_RNA_mitoBI)
seur21_RNA_mitoBI <- RunPCA(object = seur21_RNA_mitoBI, npcs = 50, verbose = FALSE)
seur21_RNA_mitoBI <- RunUMAP(object = seur21_RNA_mitoBI, reduction = "pca", dims = 1:50)
seur21_RNA_mitoBI <- FindNeighbors(seur21_RNA_mitoBI,reduction="pca",dims=1:50)
seur21_RNA_mitoBI <- FindClusters(seur21_RNA_mitoBI,resolution=0.8)

## pK Identification (no ground-truth)
sweep.res.list_seur21_RNA_mitoBI <- paramSweep_v3(seur21_RNA_mitoBI, PCs = 1:50, sct = FALSE)
sweep.stats_seur21_RNA_mitoBI <- summarizeSweep(sweep.res.list_seur21_RNA_mitoBI, GT = FALSE)
bcmvn_seur21_RNA_mitoBI <- find.pK(sweep.stats_seur21_RNA_mitoBI)
bcmvn_seur21_RNA_mitoBI

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seur21_RNA_mitoBI@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.05*nrow(seur21_RNA_mitoBI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seur21_RNA_mitoBI <- doubletFinder_v3(seur21_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.21, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
#seur21_RNA_mitoBI <- doubletFinder_v3(seur21_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
DimPlot(seur21_RNA_mitoBI, pt.size=0.1, group.by="DF.classifications_0.25_0.21_443")

## paste column to new place and delete old one
seur21_RNA_mitoBI$sample_doublets <- seur21_RNA_mitoBI$DF.classifications_0.25_0.21_443
seur21_RNA_mitoBI$DF.classifications_0.25_0.21_443 <- NULL

##remove doublets and shrink object
seur21_RNA_singlet<-subset(x=seur21_RNA_mitoBI, subset = DF.classifications_0.25_0.21_443 == "Singlet")
seur21_RNA_singlet_diet<-DietSeurat(seur21_RNA_singlet)
seur21_RNA_singlet_diet$DF.classifications_0.25_0.21_443<-NULL
seur21_RNA_singlet_diet$pANN_0.25_0.21_443<-NULL
seur21_RNA_singlet_diet$seurat_clusters<-NULL
seur21_RNA_singlet_diet$RNA_snn_res.0.8<-NULL
saveRDS(seur21_RNA_singlet_diet,file="objects/seur21_RNA.RDS")

DefaultAssay(object = seur22_RNA_mitoBI) <- "RNA"
seur22_RNA_mitoBI <- ScaleData(object = seur22_RNA_mitoBI, verbose = FALSE)
seur22_RNA_mitoBI <- FindVariableFeatures(seur22_RNA_mitoBI)
seur22_RNA_mitoBI <- RunPCA(object = seur22_RNA_mitoBI, npcs = 50, verbose = FALSE)
seur22_RNA_mitoBI <- RunUMAP(object = seur22_RNA_mitoBI, reduction = "pca", dims = 1:50)
seur22_RNA_mitoBI <- FindNeighbors(seur22_RNA_mitoBI,reduction="pca",dims=1:50)
seur22_RNA_mitoBI <- FindClusters(seur22_RNA_mitoBI,resolution=0.8)

## pK Identification (no ground-truth)
sweep.res.list_seur22_RNA_mitoBI <- paramSweep_v3(seur22_RNA_mitoBI, PCs = 1:50, sct = FALSE)
sweep.stats_seur22_RNA_mitoBI <- summarizeSweep(sweep.res.list_seur22_RNA_mitoBI, GT = FALSE)
bcmvn_seur22_RNA_mitoBI <- find.pK(sweep.stats_seur22_RNA_mitoBI)
bcmvn_seur22_RNA_mitoBI

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seur22_RNA_mitoBI@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.05*nrow(seur22_RNA_mitoBI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seur22_RNA_mitoBI <- doubletFinder_v3(seur22_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.14, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
#seur22_RNA_mitoBI <- doubletFinder_v3(seur22_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
DimPlot(seur22_RNA_mitoBI, pt.size=0.1, group.by="DF.classifications_0.25_0.14_165")

## paste column to new place and delete old one
seur22_RNA_mitoBI$sample_doublets <- seur22_RNA_mitoBI$DF.classifications_0.25_0.14_165
seur22_RNA_mitoBI$DF.classifications_0.25_0.14_165 <- NULL

##remove doublets and shrink object
seur22_RNA_singlet<-subset(x=seur22_RNA_mitoBI, subset = DF.classifications_0.25_0.14_165 == "Singlet")
seur22_RNA_singlet_diet<-DietSeurat(seur22_RNA_singlet)
seur22_RNA_singlet_diet$DF.classifications_0.25_0.14_165<-NULL
seur22_RNA_singlet_diet$pANN_0.25_0.14_165<-NULL
seur22_RNA_singlet_diet$seurat_clusters<-NULL
seur22_RNA_singlet_diet$RNA_snn_res.0.8<-NULL
saveRDS(seur22_RNA_singlet_diet,file="objects/seur22_RNA.RDS")

DefaultAssay(object = seur23_RNA_mitoBI) <- "RNA"
seur23_RNA_mitoBI <- ScaleData(object = seur23_RNA_mitoBI, verbose = FALSE)
seur23_RNA_mitoBI <- FindVariableFeatures(seur23_RNA_mitoBI)
seur23_RNA_mitoBI <- RunPCA(object = seur23_RNA_mitoBI, npcs = 50, verbose = FALSE)
seur23_RNA_mitoBI <- RunUMAP(object = seur23_RNA_mitoBI, reduction = "pca", dims = 1:50)
seur23_RNA_mitoBI <- FindNeighbors(seur23_RNA_mitoBI,reduction="pca",dims=1:50)
seur23_RNA_mitoBI <- FindClusters(seur23_RNA_mitoBI,resolution=0.8)

## pK Identification (no ground-truth)
sweep.res.list_seur23_RNA_mitoBI <- paramSweep_v3(seur23_RNA_mitoBI, PCs = 1:50, sct = FALSE)
sweep.stats_seur23_RNA_mitoBI <- summarizeSweep(sweep.res.list_seur23_RNA_mitoBI, GT = FALSE)
bcmvn_seur23_RNA_mitoBI <- find.pK(sweep.stats_seur23_RNA_mitoBI)
bcmvn_seur23_RNA_mitoBI

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seur23_RNA_mitoBI@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.05*nrow(seur23_RNA_mitoBI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seur23_RNA_mitoBI <- doubletFinder_v3(seur23_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.17, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
#seur23_RNA_mitoBI <- doubletFinder_v3(seur23_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
DimPlot(seur23_RNA_mitoBI, pt.size=0.1, group.by="DF.classifications_0.25_0.17_310")

## paste column to new place and delete old one
seur23_RNA_mitoBI$sample_doublets <- seur23_RNA_mitoBI$DF.classifications_0.25_0.17_310
seur23_RNA_mitoBI$DF.classifications_0.25_0.17_310 <- NULL

##remove doublets and shrink object
seur23_RNA_singlet<-subset(x=seur23_RNA_mitoBI, subset = DF.classifications_0.25_0.17_310 == "Singlet")
seur23_RNA_singlet_diet<-DietSeurat(seur23_RNA_singlet)
seur23_RNA_singlet_diet$DF.classifications_0.25_0.17_310<-NULL
seur23_RNA_singlet_diet$pANN_0.25_0.17_310<-NULL
seur23_RNA_singlet_diet$seurat_clusters<-NULL
seur23_RNA_singlet_diet$RNA_snn_res.0.8<-NULL
saveRDS(seur23_RNA_singlet_diet,file="objects/seur23_RNA.RDS")