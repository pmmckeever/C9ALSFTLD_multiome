library(Seurat)
library(ggplot2)
library(viridis)
library(scClustViz)
library(org.Hs.eg.db)
library(DoubletFinder)
library(future)
set
plan("multisession", workers = 8)
options(future.globals.maxSize = 120 * 1024^3)

#C9ALSFTLD matrix location
input_from_10x11 <- "~/projects/def-rogaeva/prime235/BAMfree/snRNA/RNA_02-fctx/filtered_feature_bc_matrix"
input_from_10x12 <- "~/projects/def-rogaeva/prime235/BAMfree/snRNA/RNA_14a-fctx/filtered_feature_bc_matrix"
input_from_10x13 <- "~/projects/def-rogaeva/prime235/BAMfree/snRNA/RNA_22-fctx/filtered_feature_bc_matrix"
input_from_10x14 <- "~/projects/def-rogaeva/prime235/BAMfree/snRNA/RNA_37-fctx/filtered_feature_bc_matrix"
input_from_10x15 <- "~/projects/def-rogaeva/prime235/BAMfree/snRNA/RNA_46-fctx/filtered_feature_bc_matrix"
input_from_10x16 <- "~/projects/def-rogaeva/prime235/BAMfree/snRNA/RNA_47-fctx/filtered_feature_bc_matrix"
input_from_10x17 <- "~/projects/def-rogaeva/prime235/BAMfree/snRNA/RNA_76a-fctx/filtered_feature_bc_matrix"


#C9ALSFTLD cases read-in
seur11_RNA <- CreateSeuratObject(counts=Read10X(input_from_10x11), project="C9ALSFTLD1",
                          min.cells=1,min.features=1)
seur12_RNA <- CreateSeuratObject(counts=Read10X(input_from_10x12), project="C9ALSFTLD2",
                          min.cells=1,min.features=1)
seur13_RNA <- CreateSeuratObject(counts=Read10X(input_from_10x13), project="C9ALSFTLD3",
                            min.cells=1,min.features=1)
seur14_RNA <- CreateSeuratObject(counts=Read10X(input_from_10x14), project="C9ALSFTLD4",
                            min.cells=1,min.features=1)
seur15_RNA <- CreateSeuratObject(counts=Read10X(input_from_10x15), project="C9ALSFTLD5",
                          min.cells=1,min.features=1)
seur16_RNA <- CreateSeuratObject(counts=Read10X(input_from_10x16), project="C9ALSFTLD6",
                            min.cells=1,min.features=1)
seur17_RNA <- CreateSeuratObject(counts=Read10X(input_from_10x17), project="C9ALSFTLD7",
                          min.cells=1,min.features=1)

show(seur11_RNA)
show(seur12_RNA)
show(seur13_RNA)
show(seur14_RNA)
show(seur15_RNA)
show(seur16_RNA)
show(seur17_RNA)

seur11_RNA[["percent.mt"]] <- PercentageFeatureSet(seur11_RNA, pattern = "^MT-")
seur12_RNA[["percent.mt"]] <- PercentageFeatureSet(seur12_RNA, pattern = "^MT-")
seur13_RNA[["percent.mt"]] <- PercentageFeatureSet(seur13_RNA, pattern = "^MT-")
seur14_RNA[["percent.mt"]] <- PercentageFeatureSet(seur14_RNA, pattern = "^MT-")
seur15_RNA[["percent.mt"]] <- PercentageFeatureSet(seur15_RNA, pattern = "^MT-")
seur16_RNA[["percent.mt"]] <- PercentageFeatureSet(seur16_RNA, pattern = "^MT-")
seur17_RNA[["percent.mt"]] <- PercentageFeatureSet(seur17_RNA, pattern = "^MT-")

seur11_RNA[["percent.rps"]] <- PercentageFeatureSet(seur11_RNA, pattern = "^RPS")
seur12_RNA[["percent.rps"]] <- PercentageFeatureSet(seur12_RNA, pattern = "^RPS")
seur13_RNA[["percent.rps"]] <- PercentageFeatureSet(seur13_RNA, pattern = "^RPS")
seur14_RNA[["percent.rps"]] <- PercentageFeatureSet(seur14_RNA, pattern = "^RPS")
seur15_RNA[["percent.rps"]] <- PercentageFeatureSet(seur15_RNA, pattern = "^RPS")
seur16_RNA[["percent.rps"]] <- PercentageFeatureSet(seur16_RNA, pattern = "^RPS")
seur17_RNA[["percent.rps"]] <- PercentageFeatureSet(seur17_RNA, pattern = "^RPS")

seur11_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur11_RNA, pattern = "^RPL")
seur12_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur12_RNA, pattern = "^RPL")
seur13_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur13_RNA, pattern = "^RPL")
seur14_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur14_RNA, pattern = "^RPL")
seur15_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur15_RNA, pattern = "^RPL")
seur16_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur16_RNA, pattern = "^RPL")
seur17_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur17_RNA, pattern = "^RPL")

mito_gene_identifier <- "^MT-"
mads_thresh <- 3
hard_thresh <- 50
seur11_RNA <- PercentageFeatureSet(seur11_RNA, pattern = "^MT-", col.name = "percent.mt")
mito_thresh <- median(seur11_RNA$percent.mt) + mad(seur11_RNA$percent.mt) * mads_thresh
drop_mito <- seur11_RNA$percent.mt > mito_thresh | seur11_RNA$percent.mt > hard_thresh
seur11_RNA_mitoBI <- seur11_RNA[,!drop_mito]
show(seur11_RNA_mitoBI)
seur12_RNA <- PercentageFeatureSet(seur12_RNA, pattern = "^MT-", col.name = "percent.mt")
mito_thresh <- median(seur12_RNA$percent.mt) + mad(seur12_RNA$percent.mt) * mads_thresh
drop_mito <- seur12_RNA$percent.mt > mito_thresh | seur12_RNA$percent.mt > hard_thresh
seur12_RNA_mitoBI <- seur12_RNA[,!drop_mito]
show(seur12_RNA_mitoBI)
seur13_RNA <- PercentageFeatureSet(seur13_RNA, pattern = "^MT-", col.name = "percent.mt")
mito_thresh <- median(seur13_RNA$percent.mt) + mad(seur13_RNA$percent.mt) * mads_thresh
drop_mito <- seur13_RNA$percent.mt > mito_thresh | seur13_RNA$percent.mt > hard_thresh
seur13_RNA_mitoBI <- seur13_RNA[,!drop_mito]
show(seur13_RNA_mitoBI)
seur14_RNA <- PercentageFeatureSet(seur14_RNA, pattern = "^MT-", col.name = "percent.mt")
mito_thresh <- median(seur14_RNA$percent.mt) + mad(seur14_RNA$percent.mt) * mads_thresh
drop_mito <- seur14_RNA$percent.mt > mito_thresh | seur14_RNA$percent.mt > hard_thresh
seur14_RNA_mitoBI <- seur14_RNA[,!drop_mito]
show(seur14_RNA_mitoBI)
seur15_RNA <- PercentageFeatureSet(seur15_RNA, pattern = "^MT-", col.name = "percent.mt")
mito_thresh <- median(seur15_RNA$percent.mt) + mad(seur15_RNA$percent.mt) * mads_thresh
drop_mito <- seur15_RNA$percent.mt > mito_thresh | seur15_RNA$percent.mt > hard_thresh
seur15_RNA_mitoBI <- seur15_RNA[,!drop_mito]
show(seur15_RNA_mitoBI)
seur16_RNA <- PercentageFeatureSet(seur16_RNA, pattern = "^MT-", col.name = "percent.mt")
mito_thresh <- median(seur16_RNA$percent.mt) + mad(seur16_RNA$percent.mt) * mads_thresh
drop_mito <- seur16_RNA$percent.mt > mito_thresh | seur16_RNA$percent.mt > hard_thresh
seur16_RNA_mitoBI <- seur16_RNA[,!drop_mito]
show(seur16_RNA_mitoBI)
seur17_RNA <- PercentageFeatureSet(seur17_RNA, pattern = "^MT-", col.name = "percent.mt")
mito_thresh <- median(seur17_RNA$percent.mt) + mad(seur17_RNA$percent.mt) * mads_thresh
drop_mito <- seur17_RNA$percent.mt > mito_thresh | seur17_RNA$percent.mt > hard_thresh
seur17_RNA_mitoBI <- seur17_RNA[,!drop_mito]
show(seur17_RNA_mitoBI)

seur11_RNA_mitoBI <- subset(seur11_RNA_mitoBI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)
seur12_RNA_mitoBI <- subset(seur12_RNA_mitoBI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)
seur13_RNA_mitoBI <- subset(seur13_RNA_mitoBI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)
seur14_RNA_mitoBI <- subset(seur14_RNA_mitoBI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)
seur15_RNA_mitoBI <- subset(seur15_RNA_mitoBI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)
seur16_RNA_mitoBI <- subset(seur16_RNA_mitoBI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)
seur17_RNA_mitoBI <- subset(seur17_RNA_mitoBI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)

show(seur11_RNA_mitoBI)
show(seur12_RNA_mitoBI)    
show(seur13_RNA_mitoBI)
show(seur14_RNA_mitoBI)
show(seur15_RNA_mitoBI)
show(seur16_RNA_mitoBI)
show(seur17_RNA_mitoBI)

seur11_RNA_mitoBI$sample <- 'C9ALSFTLD1'
seur12_RNA_mitoBI$sample <- 'C9ALSFTLD2'
seur13_RNA_mitoBI$sample <- 'C9ALSFTLD3'
seur14_RNA_mitoBI$sample <- 'C9ALSFTLD4'
seur15_RNA_mitoBI$sample <- 'C9ALSFTLD5'
seur16_RNA_mitoBI$sample <- 'C9ALSFTLD6'
seur17_RNA_mitoBI$sample <- 'C9ALSFTLD7'

seur11_RNA_mitoBI<-NormalizeData(seur11_RNA_mitoBI)
seur12_RNA_mitoBI<-NormalizeData(seur12_RNA_mitoBI)    
seur13_RNA_mitoBI<-NormalizeData(seur13_RNA_mitoBI)
seur14_RNA_mitoBI<-NormalizeData(seur14_RNA_mitoBI)
seur15_RNA_mitoBI<-NormalizeData(seur15_RNA_mitoBI)
seur16_RNA_mitoBI<-NormalizeData(seur16_RNA_mitoBI)
seur17_RNA_mitoBI<-NormalizeData(seur17_RNA_mitoBI)

seur11_RNA_mitoBI<-FindVariableFeatures(seur11_RNA_mitoBI)
seur12_RNA_mitoBI<-FindVariableFeatures(seur12_RNA_mitoBI)
seur13_RNA_mitoBI<-FindVariableFeatures(seur13_RNA_mitoBI)
seur14_RNA_mitoBI<-FindVariableFeatures(seur14_RNA_mitoBI)
seur15_RNA_mitoBI<-FindVariableFeatures(seur15_RNA_mitoBI)
seur16_RNA_mitoBI<-FindVariableFeatures(seur16_RNA_mitoBI)
seur17_RNA_mitoBI<-FindVariableFeatures(seur17_RNA_mitoBI)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

seur11_RNA_mitoBI <- CellCycleScoring(seur11_RNA_mitoBI, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seur12_RNA_mitoBI <- CellCycleScoring(seur12_RNA_mitoBI, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seur13_RNA_mitoBI <- CellCycleScoring(seur13_RNA_mitoBI, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seur14_RNA_mitoBI <- CellCycleScoring(seur14_RNA_mitoBI, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seur15_RNA_mitoBI <- CellCycleScoring(seur15_RNA_mitoBI, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seur16_RNA_mitoBI <- CellCycleScoring(seur16_RNA_mitoBI, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seur17_RNA_mitoBI <- CellCycleScoring(seur17_RNA_mitoBI, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DefaultAssay(object = seur11_RNA_mitoBI) <- "RNA"
seur11_RNA_mitoBI <- ScaleData(object = seur11_RNA_mitoBI, verbose = FALSE)
seur11_RNA_mitoBI <- FindVariableFeatures(seur11_RNA_mitoBI)
seur11_RNA_mitoBI <- RunPCA(object = seur11_RNA_mitoBI, npcs = 50, verbose = FALSE)
seur11_RNA_mitoBI <- RunUMAP(object = seur11_RNA_mitoBI, reduction = "pca", dims = 1:50)
seur11_RNA_mitoBI <- FindNeighbors(seur11_RNA_mitoBI,reduction="pca",dims=1:50)
seur11_RNA_mitoBI <- FindClusters(seur11_RNA_mitoBI,resolution=0.8)

## pK Identification (no ground-truth)
sweep.res.list_seur11_RNA_mitoBI <- paramSweep_v3(seur11_RNA_mitoBI, PCs = 1:50, sct = FALSE)
sweep.stats_seur11_RNA_mitoBI <- summarizeSweep(sweep.res.list_seur11_RNA_mitoBI, GT = FALSE)
bcmvn_seur11_RNA_mitoBI <- find.pK(sweep.stats_seur11_RNA_mitoBI)
bcmvn_seur11_RNA_mitoBI

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seur11_RNA_mitoBI@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.05*nrow(seur11_RNA_mitoBI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seur11_RNA_mitoBI <- doubletFinder_v3(seur11_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.3, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
#seur11_RNA_mitoBI <- doubletFinder_v3(seur11_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
DimPlot(seur11_RNA_mitoBI, pt.size=0.1, group.by="DF.classifications_0.25_0.3_202")

## paste column to new place and delete old one
seur11_RNA_mitoBI$sample_doublets <- seur11_RNA_mitoBI$DF.classifications_0.25_0.3_202
seur11_RNA_mitoBI$DF.classifications_0.25_0.3_202 <- NULL

##remove doublets and shrink object
seur11_RNA_singlet<-subset(x=seur11_RNA_mitoBI, subset = DF.classifications_0.25_0.3_202 == "Singlet")
seur11_RNA_singlet_diet<-DietSeurat(seur11_RNA_singlet)
seur11_RNA_singlet_diet$DF.classifications_0.25_0.3_202<-NULL
seur11_RNA_singlet_diet$pANN_0.25_0.3_202<-NULL
seur11_RNA_singlet_diet$seurat_clusters<-NULL
seur11_RNA_singlet_diet$RNA_snn_res.0.8<-NULL
saveRDS(seur11_RNA_singlet_diet,file="objects/seur11_RNA.RDS")

DefaultAssay(object = seur12_RNA_mitoBI) <- "RNA"
seur12_RNA_mitoBI <- ScaleData(object = seur12_RNA_mitoBI, verbose = FALSE)
seur12_RNA_mitoBI <- FindVariableFeatures(seur12_RNA_mitoBI)
seur12_RNA_mitoBI <- RunPCA(object = seur12_RNA_mitoBI, npcs = 50, verbose = FALSE)
seur12_RNA_mitoBI <- RunUMAP(object = seur12_RNA_mitoBI, reduction = "pca", dims = 1:50)
seur12_RNA_mitoBI <- FindNeighbors(seur12_RNA_mitoBI,reduction="pca",dims=1:50)
seur12_RNA_mitoBI <- FindClusters(seur12_RNA_mitoBI,resolution=0.8)

## pK Identification (no ground-truth)
sweep.res.list_seur12_RNA_mitoBI <- paramSweep_v3(seur12_RNA_mitoBI, PCs = 1:50, sct = FALSE)
sweep.stats_seur12_RNA_mitoBI <- summarizeSweep(sweep.res.list_seur12_RNA_mitoBI, GT = FALSE)
bcmvn_seur12_RNA_mitoBI <- find.pK(sweep.stats_seur12_RNA_mitoBI)
bcmvn_seur12_RNA_mitoBI

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seur12_RNA_mitoBI@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.05*nrow(seur12_RNA_mitoBI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seur12_RNA_mitoBI <- doubletFinder_v3(seur12_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.29, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
#seur12_RNA_mitoBI <- doubletFinder_v3(seur12_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
DimPlot(seur12_RNA_mitoBI, pt.size=0.1, group.by="DF.classifications_0.25_0.29_391")

## paste column to new place and delete old one
seur12_RNA_mitoBI$sample_doublets <- seur12_RNA_mitoBI$DF.classifications_0.25_0.29_391
seur12_RNA_mitoBI$DF.classifications_0.25_0.29_391 <- NULL

##remove doublets and shrink object
seur12_RNA_singlet<-subset(x=seur12_RNA_mitoBI, subset = DF.classifications_0.25_0.29_391 == "Singlet")
seur12_RNA_singlet_diet<-DietSeurat(seur12_RNA_singlet)
seur12_RNA_singlet_diet$DF.classifications_0.25_0.29_391<-NULL
seur12_RNA_singlet_diet$pANN_0.25_0.29_391<-NULL
seur12_RNA_singlet_diet$seurat_clusters<-NULL
seur12_RNA_singlet_diet$RNA_snn_res.0.8<-NULL
saveRDS(seur12_RNA_singlet_diet,file="objects/seur12_RNA.RDS")

DefaultAssay(object = seur13_RNA_mitoBI) <- "RNA"
seur13_RNA_mitoBI <- ScaleData(object = seur13_RNA_mitoBI, verbose = FALSE)
seur13_RNA_mitoBI <- FindVariableFeatures(seur13_RNA_mitoBI)
seur13_RNA_mitoBI <- RunPCA(object = seur13_RNA_mitoBI, npcs = 50, verbose = FALSE)
seur13_RNA_mitoBI <- RunUMAP(object = seur13_RNA_mitoBI, reduction = "pca", dims = 1:50)
seur13_RNA_mitoBI <- FindNeighbors(seur13_RNA_mitoBI,reduction="pca",dims=1:50)
seur13_RNA_mitoBI <- FindClusters(seur13_RNA_mitoBI,resolution=0.8)

## pK Identification (no ground-truth)
sweep.res.list_seur13_RNA_mitoBI <- paramSweep_v3(seur13_RNA_mitoBI, PCs = 1:50, sct = FALSE)
sweep.stats_seur13_RNA_mitoBI <- summarizeSweep(sweep.res.list_seur13_RNA_mitoBI, GT = FALSE)
bcmvn_seur13_RNA_mitoBI <- find.pK(sweep.stats_seur13_RNA_mitoBI)
bcmvn_seur13_RNA_mitoBI

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seur13_RNA_mitoBI@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.05*nrow(seur13_RNA_mitoBI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seur13_RNA_mitoBI <- doubletFinder_v3(seur13_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.18, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
#seur13_RNA_mitoBI <- doubletFinder_v3(seur13_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
DimPlot(seur13_RNA_mitoBI, pt.size=0.1, group.by="DF.classifications_0.25_0.18_405")

## paste column to new place and delete old one
seur13_RNA_mitoBI$sample_doublets <- seur13_RNA_mitoBI$DF.classifications_0.25_0.18_405
seur13_RNA_mitoBI$DF.classifications_0.25_0.18_405 <- NULL

##remove doublets and shrink object
seur13_RNA_singlet<-subset(x=seur13_RNA_mitoBI, subset = DF.classifications_0.25_0.18_405 == "Singlet")
seur13_RNA_singlet_diet<-DietSeurat(seur13_RNA_singlet)
seur13_RNA_singlet_diet$DF.classifications_0.25_0.18_405<-NULL
seur13_RNA_singlet_diet$pANN_0.25_0.18_405<-NULL
seur13_RNA_singlet_diet$seurat_clusters<-NULL
seur13_RNA_singlet_diet$RNA_snn_res.0.8<-NULL
saveRDS(seur13_RNA_singlet_diet,file="objects/seur13_RNA.RDS")

DefaultAssay(object = seur14_RNA_mitoBI) <- "RNA"
seur14_RNA_mitoBI <- ScaleData(object = seur14_RNA_mitoBI, verbose = FALSE)
seur14_RNA_mitoBI <- FindVariableFeatures(seur14_RNA_mitoBI)
seur14_RNA_mitoBI <- RunPCA(object = seur14_RNA_mitoBI, npcs = 50, verbose = FALSE)
seur14_RNA_mitoBI <- RunUMAP(object = seur14_RNA_mitoBI, reduction = "pca", dims = 1:50)
seur14_RNA_mitoBI <- FindNeighbors(seur14_RNA_mitoBI,reduction="pca",dims=1:50)
seur14_RNA_mitoBI <- FindClusters(seur14_RNA_mitoBI,resolution=0.8)

# pK Identification (no ground-truth)
sweep.res.list_seur14_RNA_mitoBI <- paramSweep_v3(seur14_RNA_mitoBI, PCs = 1:50, sct = FALSE)
sweep.stats_seur14_RNA_mitoBI <- summarizeSweep(sweep.res.list_seur14_RNA_mitoBI, GT = FALSE)
bcmvn_seur14_RNA_mitoBI <- find.pK(sweep.stats_seur14_RNA_mitoBI)
bcmvn_seur14_RNA_mitoBI

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seur14_RNA_mitoBI@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.05*nrow(seur14_RNA_mitoBI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seur14_RNA_mitoBI <- doubletFinder_v3(seur14_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.24, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
#seur14_RNA_mitoBI <- doC9_ALS_FTLD_RNA,file="Jan21_C9ALS_FTLD_RNA.RDS"ubletFinder_v3(seur14_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
DimPlot(seur14_RNA_mitoBI, pt.size=0.1, group.by="DF.classifications_0.25_0.26_225")

## paste column to new place and delete old one
seur14_RNA_mitoBI$sample_doublets <- seur14_RNA_mitoBI$DF.classifications_0.25_0.24_233
seur14_RNA_mitoBI$DF.classifications_0.25_0.24_233 <- NULL

##remove doublets and shrink object
seur14_RNA_singlet<-subset(x=seur14_RNA_mitoBI, subset = DF.classifications_0.25_0.24_233 == "Singlet")
seur14_RNA_singlet_diet<-DietSeurat(seur14_RNA_singlet)
seur14_RNA_singlet_diet$DF.classifications_0.25_0.24_233<-NULL
seur14_RNA_singlet_diet$pANN_0.25_0.24_233<-NULL
seur14_RNA_singlet_diet$seurat_clusters<-NULL
seur14_RNA_singlet_diet$RNA_snn_res.0.8<-NULL
saveRDS(seur14_RNA_singlet_diet,file="objects/seur14_RNA.RDS")

DefaultAssay(object = seur15_RNA_mitoBI) <- "RNA"
seur15_RNA_mitoBI <- ScaleData(object = seur15_RNA_mitoBI, verbose = FALSE)
seur15_RNA_mitoBI <- FindVariableFeatures(seur15_RNA_mitoBI)
seur15_RNA_mitoBI <- RunPCA(object = seur15_RNA_mitoBI, npcs = 50, verbose = FALSE)
seur15_RNA_mitoBI <- RunUMAP(object = seur15_RNA_mitoBI, reduction = "pca", dims = 1:50)
seur15_RNA_mitoBI <- FindNeighbors(seur15_RNA_mitoBI,reduction="pca",dims=1:50)
seur15_RNA_mitoBI <- FindClusters(seur15_RNA_mitoBI,resolution=0.8)

# pK Identification (no ground-truth)
sweep.res.list_seur15_RNA_mitoBI <- paramSweep_v3(seur15_RNA_mitoBI, PCs = 1:50, sct = FALSE)
sweep.stats_seur15_RNA_mitoBI <- summarizeSweep(sweep.res.list_seur15_RNA_mitoBI, GT = FALSE)
bcmvn_seur15_RNA_mitoBI <- find.pK(sweep.stats_seur15_RNA_mitoBI)
bcmvn_seur15_RNA_mitoBI

# Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seur15_RNA_mitoBI@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.05*nrow(seur15_RNA_mitoBI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seur15_RNA_mitoBI <- doubletFinder_v3(seur15_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.23, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
#seur15_RNA_mitoBI <- doubletFinder_v3(seur15_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
DimPlot(seur15_RNA_mitoBI, pt.size=0.1, group.by="DF.classifications_0.25_0.21_184")

## paste column to new place and delete old one
seur15_RNA_mitoBI$sample_doublets <- seur15_RNA_mitoBI$DF.classifications_0.25_0.23_164
seur15_RNA_mitoBI$DF.classifications_0.25_0.23_164 <- NULL

##remove doublets and shrink object
seur15_RNA_singlet<-subset(x=seur15_RNA_mitoBI, subset = DF.classifications_0.25_0.23_164 == "Singlet")
seur15_RNA_singlet_diet<-DietSeurat(seur15_RNA_singlet)
seur15_RNA_singlet_diet$DF.classifications_0.25_0.23_164<-NULL
seur15_RNA_singlet_diet$pANN_0.25_0.23_164<-NULL
seur15_RNA_singlet_diet$seurat_clusters<-NULL
seur15_RNA_singlet_diet$RNA_snn_res.0.8<-NULL
saveRDS(seur15_RNA_singlet_diet,file="objects/seur15_RNA.RDS")

DefaultAssay(object = seur16_RNA_mitoBI) <- "RNA"
seur16_RNA_mitoBI <- ScaleData(object = seur16_RNA_mitoBI, verbose = FALSE)
seur16_RNA_mitoBI <- FindVariableFeatures(seur16_RNA_mitoBI)
seur16_RNA_mitoBI <- RunPCA(object = seur16_RNA_mitoBI, npcs = 50, verbose = FALSE)
seur16_RNA_mitoBI <- RunUMAP(object = seur16_RNA_mitoBI, reduction = "pca", dims = 1:50)
seur16_RNA_mitoBI <- FindNeighbors(seur16_RNA_mitoBI,reduction="pca",dims=1:50)
seur16_RNA_mitoBI <- FindClusters(seur16_RNA_mitoBI,resolution=0.8)

# pK Identification (no ground-truth)
sweep.res.list_seur16_RNA_mitoBI <- paramSweep_v3(seur16_RNA_mitoBI, PCs = 1:50, sct = FALSE)
sweep.stats_seur16_RNA_mitoBI <- summarizeSweep(sweep.res.list_seur16_RNA_mitoBI, GT = FALSE)
bcmvn_seur16_RNA_mitoBI <- find.pK(sweep.stats_seur16_RNA_mitoBI)
bcmvn_seur16_RNA_mitoBI

# Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seur16_RNA_mitoBI@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.05*nrow(seur16_RNA_mitoBI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seur16_RNA_mitoBI <- doubletFinder_v3(seur16_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
#seur16_RNA_mitoBI <- doubletFinder_v3(seur16_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
DimPlot(seur16_RNA_mitoBI, pt.size=0.1, group.by="DF.classifications_0.25_0.01_188")

## paste column to new place and delete old one
seur16_RNA_mitoBI$sample_doublets <- seur16_RNA_mitoBI$DF.classifications_0.25_0.01_188
seur16_RNA_mitoBI$DF.classifications_0.25_0.01_188 <- NULL

##remove doublets and shrink object
seur16_RNA_singlet<-subset(x=seur16_RNA_mitoBI, subset = DF.classifications_0.25_0.01_188 == "Singlet")
seur16_RNA_singlet_diet<-DietSeurat(seur16_RNA_singlet)
seur16_RNA_singlet_diet$DF.classifications_0.25_0.01_188<-NULL
seur16_RNA_singlet_diet$pANN_0.25_0.01_188<-NULL
seur16_RNA_singlet_diet$seurat_clusters<-NULL
seur16_RNA_singlet_diet$RNA_snn_res.0.8<-NULL
saveRDS(seur16_RNA_singlet_diet,file="objects/seur16_RNA.RDS")

DefaultAssay(object = seur17_RNA_mitoBI) <- "RNA"
seur17_RNA_mitoBI <- ScaleData(object = seur17_RNA_mitoBI, verbose = FALSE)
seur17_RNA_mitoBI <- FindVariableFeatures(seur17_RNA_mitoBI)
seur17_RNA_mitoBI <- RunPCA(object = seur17_RNA_mitoBI, npcs = 50, verbose = FALSE)
seur17_RNA_mitoBI <- RunUMAP(object = seur17_RNA_mitoBI, reduction = "pca", dims = 1:50)
seur17_RNA_mitoBI <- FindNeighbors(seur17_RNA_mitoBI,reduction="pca",dims=1:50)
seur17_RNA_mitoBI <- FindClusters(seur17_RNA_mitoBI,resolution=0.8)

## pK Identification (no ground-truth)
sweep.res.list_seur17_RNA_mitoBI <- paramSweep_v3(seur17_RNA_mitoBI, PCs = 1:50, sct = FALSE)
sweep.stats_seur17_RNA_mitoBI <- summarizeSweep(sweep.res.list_seur17_RNA_mitoBI, GT = FALSE)
bcmvn_seur17_RNA_mitoBI <- find.pK(sweep.stats_seur17_RNA_mitoBI)

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seur17_RNA_mitoBI@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)         
nExp_poi <- round(0.05*nrow(seur17_RNA_mitoBI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seur17_RNA_mitoBI <- doubletFinder_v3(seur17_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
#seur17_RNA_mitoBI <- doubletFinder_v3(seur17_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
DimPlot(seur17_RNA_mitoBI, pt.size=0.1, group.by="DF.classifications_0.25_0.005_341")

## paste column to new place and delete old one
seur17_RNA_mitoBI$sample_doublets <- seur17_RNA_mitoBI$DF.classifications_0.25_0.005_341
seur17_RNA_mitoBI$DF.classifications_0.25_0.005_341 <- NULL

##remove doublets and shrink object
seur17_RNA_singlet<-subset(x=seur17_RNA_mitoBI, subset = DF.classifications_0.25_0.005_341== "Singlet")
seur17_RNA_singlet_diet<-DietSeurat(seur17_RNA_singlet)
seur17_RNA_singlet_diet$DF.classifications_0.25_0.005_341<-NULL
seur17_RNA_singlet_diet$pANN_0.25_0.005_341<-NULL
seur17_RNA_singlet_diet$seurat_clusters<-NULL
seur17_RNA_singlet_diet$RNA_snn_res.0.8<-NULL
saveRDS(seur17_RNA_singlet_diet,file="objects/seur17_RNA.RDS")