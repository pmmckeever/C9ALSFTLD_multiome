library(Seurat)
heimlibrary(org.Hs.eg.db)
library(DoubletFinder)
set.seed(1234)

#CTRL matrix location
input_from_10xC1 <- "/ALS_snRNA/CTRL1_snRNA/outs/filtered_feature_bc_matrix"
input_from_10xC2 <- "/ALS_snRNA/CTRL2_snRNA/outs/filtered_feature_bc_matrix"
input_from_10xC3 <- "/ALS_snRNA/CTRL3_snRNA/outs/filtered_feature_bc_matrix"
input_from_10xC4 <- "/ALS_snRNA/CTRL4_snRNA/outs/filtered_feature_bc_matrix"
input_from_10xC5 <- "/ALS_snRNA/CTRL5_snRNA/outs/filtered_feature_bc_matrix"
input_from_10xC6 <- "/ALS_snRNA/CTRL6_snRNA/outs/filtered_feature_bc_matrix"

#control cases read-in
seurC1_RNA <- CreateSeuratObject(counts=Read10X(input_from_10xC1), project="CTRL1",
                          min.cells=1,min.features=1)
seurC2_RNA <- CreateSeuratObject(counts=Read10X(input_from_10xC2), project="CTRL2",
                          min.cells=1,min.features=1)
seurC3_RNA <- CreateSeuratObject(counts=Read10X(input_from_10xC3), project="CTRL3",
                            min.cells=1,min.features=1)
seurC4_RNA <- CreateSeuratObject(counts=Read10X(input_from_10xC4), project="CTRL4",
                            min.cells=1,min.features=1)
seurC5_RNA <- CreateSeuratObject(counts=Read10X(input_from_10xC5), project="CTRL5",
                          min.cells=1,min.features=1)
seurC6_RNA <- CreateSeuratObject(counts=Read10X(input_from_10xC6), project="CTRL6",
                            min.cells=1,min.features=1)

show(seurC1_RNA)
show(seurC2_RNA)
show(seurC3_RNA)
show(seurC4_RNA)
show(seurC5_RNA)
show(seurC6_RNA)

seurC1_RNA[["percent.mt"]] <- PercentageFeatureSet(seurC1_RNA, pattern = "^MT-")
seurC2_RNA[["percent.mt"]] <- PercentageFeatureSet(seurC2_RNA, pattern = "^MT-")
seurC3_RNA[["percent.mt"]] <- PercentageFeatureSet(seurC3_RNA, pattern = "^MT-")
seurC4_RNA[["percent.mt"]] <- PercentageFeatureSet(seurC4_RNA, pattern = "^MT-")
seurC5_RNA[["percent.mt"]] <- PercentageFeatureSet(seurC5_RNA, pattern = "^MT-")
seurC6_RNA[["percent.mt"]] <- PercentageFeatureSet(seurC6_RNA, pattern = "^MT-")

mito_gene_identifier <- "^MT-"
mads_thresh <- 3
hard_thresh <- 50
seurC1_RNA <- PercentageFeatureSet(seurC1_RNA, pattern = "^MT-", col.name = "percent.mt")
mito_thresh <- median(seurC1_RNA$percent.mt) + mad(seurC1_RNA$percent.mt) * mads_thresh
drop_mito <- seurC1_RNA$percent.mt > mito_thresh | seurC1_RNA$percent.mt > hard_thresh
seurC1_RNA_mitoBI <- seurC1_RNA[,!drop_mito]
show(seurC1_RNA_mitoBI)
seurC2_RNA <- PercentageFeatureSet(seurC2_RNA, pattern = "^MT-", col.name = "percent.mt")
mito_thresh <- median(seurC2_RNA$percent.mt) + mad(seurC2_RNA$percent.mt) * mads_thresh
drop_mito <- seurC2_RNA$percent.mt > mito_thresh | seurC2_RNA$percent.mt > hard_thresh
seurC2_RNA_mitoBI <- seurC2_RNA[,!drop_mito]
show(seurC2_RNA_mitoBI)
seurC3_RNA <- PercentageFeatureSet(seurC3_RNA, pattern = "^MT-", col.name = "percent.mt")
mito_thresh <- median(seurC3_RNA$percent.mt) + mad(seurC3_RNA$percent.mt) * mads_thresh
drop_mito <- seurC3_RNA$percent.mt > mito_thresh | seurC3_RNA$percent.mt > hard_thresh
seurC3_RNA_mitoBI <- seurC3_RNA[,!drop_mito]
show(seurC3_RNA_mitoBI)
seurC4_RNA <- PercentageFeatureSet(seurC4_RNA, pattern = "^MT-", col.name = "percent.mt")
mito_thresh <- median(seurC4_RNA$percent.mt) + mad(seurC4_RNA$percent.mt) * mads_thresh
drop_mito <- seurC4_RNA$percent.mt > mito_thresh | seurC4_RNA$percent.mt > hard_thresh
seurC4_RNA_mitoBI <- seurC4_RNA[,!drop_mito]
show(seurC4_RNA_mitoBI)
seurC5_RNA <- PercentageFeatureSet(seurC5_RNA, pattern = "^MT-", col.name = "percent.mt")
mito_thresh <- median(seurC5_RNA$percent.mt) + mad(seurC5_RNA$percent.mt) * mads_thresh
drop_mito <- seurC5_RNA$percent.mt > mito_thresh | seurC5_RNA$percent.mt > hard_thresh
seurC5_RNA_mitoBI <- seurC5_RNA[,!drop_mito]
show(seurC5_RNA_mitoBI)
seurC6_RNA <- PercentageFeatureSet(seurC6_RNA, pattern = "^MT-", col.name = "percent.mt")
mito_thresh <- median(seurC6_RNA$percent.mt) + mad(seurC6_RNA$percent.mt) * mads_thresh
drop_mito <- seurC6_RNA$percent.mt > mito_thresh | seurC6_RNA$percent.mt > hard_thresh
seurC6_RNA_mitoBI <- seurC6_RNA[,!drop_mito]
show(seurC6_RNA_mitoBI)

seurC1_RNA_mitoBI <- subset(seurC1_RNA_mitoBI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)
seurC2_RNA_mitoBI <- subset(seurC2_RNA_mitoBI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)
seurC3_RNA_mitoBI <- subset(seurC3_RNA_mitoBI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)
seurC4_RNA_mitoBI <- subset(seurC4_RNA_mitoBI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)
seurC5_RNA_mitoBI <- subset(seurC5_RNA_mitoBI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)
seurC6_RNA_mitoBI <- subset(seurC6_RNA_mitoBI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)

show(seurC1_RNA_mitoBI)
show(seurC2_RNA_mitoBI)    
show(seurC3_RNA_mitoBI)
show(seurC4_RNA_mitoBI)
show(seurC5_RNA_mitoBI)
show(seurC6_RNA_mitoBI)

seurC1_RNA_mitoBI$sample <- 'CTRL1'
seurC2_RNA_mitoBI$sample <- 'CTRL2'
seurC3_RNA_mitoBI$sample <- 'CTRL3'
seurC4_RNA_mitoBI$sample <- 'CTRL4'
seurC5_RNA_mitoBI$sample <- 'CTRL5'
seurC6_RNA_mitoBI$sample <- 'CTRL6'

show(seurC1_RNA_mitoBI)
show(seurC2_RNA_mitoBI)    
show(seurC3_RNA_mitoBI)
show(seurC4_RNA_mitoBI)
show(seurC5_RNA_mitoBI)
show(seurC6_RNA_mitoBI)

seurC1_RNA_mitoBI<-NormalizeData(seurC1_RNA_mitoBI)
seurC2_RNA_mitoBI<-NormalizeData(seurC2_RNA_mitoBI)    
seurC3_RNA_mitoBI<-NormalizeData(seurC3_RNA_mitoBI)
seurC4_RNA_mitoBI<-NormalizeData(seurC4_RNA_mitoBI)
seurC5_RNA_mitoBI<-NormalizeData(seurC5_RNA_mitoBI)
seurC6_RNA_mitoBI<-NormalizeData(seurC6_RNA_mitoBI)

seurC1_RNA_mitoBI<-FindVariableFeatures(seurC1_RNA_mitoBI)
seurC2_RNA_mitoBI<-FindVariableFeatures(seurC2_RNA_mitoBI)    
seurC3_RNA_mitoBI<-FindVariableFeatures(seurC3_RNA_mitoBI)
seurC4_RNA_mitoBI<-FindVariableFeatures(seurC4_RNA_mitoBI)
seurC5_RNA_mitoBI<-FindVariableFeatures(seurC5_RNA_mitoBI)
seurC6_RNA_mitoBI<-FindVariableFeatures(seurC6_RNA_mitoBI)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

seurC1_RNA_mitoBI <- CellCycleScoring(seurC1_RNA_mitoBI, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seurC2_RNA_mitoBI <- CellCycleScoring(seurC2_RNA_mitoBI, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seurC3_RNA_mitoBI <- CellCycleScoring(seurC3_RNA_mitoBI, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seurC4_RNA_mitoBI <- CellCycleScoring(seurC4_RNA_mitoBI, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seurC5_RNA_mitoBI <- CellCycleScoring(seurC5_RNA_mitoBI, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seurC6_RNA_mitoBI <- CellCycleScoring(seurC6_RNA_mitoBI, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

### Individual sample doublets
DefaultAssay(object = seurC1_RNA_mitoBI) <- "RNA"
seurC1_RNA_mitoBI <- ScaleData(object = seurC1_RNA_mitoBI, verbose = FALSE)
seurC1_RNA_mitoBI <- FindVariableFeatures(seurC1_RNA_mitoBI)
seurC1_RNA_mitoBI <- RunPCA(object = seurC1_RNA_mitoBI, npcs = 50, verbose = FALSE)
seurC1_RNA_mitoBI <- RunUMAP(object = seurC1_RNA_mitoBI, reduction = "pca", dims = 1:50)
seurC1_RNA_mitoBI <- FindNeighbors(seurC1_RNA_mitoBI,reduction="pca",dims=1:50)
seurC1_RNA_mitoBI <- FindClusters(seurC1_RNA_mitoBI,resolution=0.8)

## pK Identification (no ground-truth)
sweep.res.list_seurC1_RNA_mitoBI <- paramSweep_v3(seurC1_RNA_mitoBI, PCs = 1:50, sct = FALSE)
sweep.stats_seurC1_RNA_mitoBI <- summarizeSweep(sweep.res.list_seurC1_RNA_mitoBI, GT = FALSE)
bcmvn_seurC1_RNA_mitoBI <- find.pK(sweep.stats_seurC1_RNA_mitoBI)
bcmvn_seurC1_RNA_mitoBI

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seurC1_RNA_mitoBI@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seurC1_RNA_mitoBI@meta.data$ClusteringResults
nExp_poi <- round(0.05*nrow(seurC1_RNA_mitoBI@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seurC1_RNA_mitoBI <- doubletFinder_v3(seurC1_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.24, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
#seurC1_RNA_mitoBI <- doubletFinder_v3(seurC1_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
DimPlot(seurC1_RNA_mitoBI, pt.size=0.1, group.by="DF.classifications_0.25_0.24_128")

## paste column to new place and delete old one
#seurC1_RNA_mitoBI$sample_doublets <- seurC1_RNA_mitoBI$DF.classifications_0.25_0.24_128
#seurC1_RNA_mitoBI$DF.classifications_0.25_0.24_128 <- NULL

##remove doublets and shrink object
seurC1_RNA_singlet<-subset(x=seurC1_RNA_mitoBI, subset = sample.doublets == "Singlet")
seurC1_RNA_singlet_diet<-DietSeurat(seurC1_RNA_singlet)
seurC1_RNA_singlet_diet$DF.classifications_0.25_0.24_128<-NULL
seurC1_RNA_singlet_diet$pANN_0.25_0.24_128<-NULL
seurC1_RNA_singlet_diet$seurat_clusters<-NULL
seurC1_RNA_singlet_diet$RNA_snn_res.0.8<-NULL
saveRDS(seurC1_RNA_singlet_diet,file="objects/seurC1_RNA.RDS")

DefaultAssay(object = seurC2_RNA_mitoBI) <- "RNA"
seurC2_RNA_mitoBI <- ScaleData(object = seurC2_RNA_mitoBI, verbose = FALSE)
seurC2_RNA_mitoBI <- FindVariableFeatures(seurC2_RNA_mitoBI)
seurC2_RNA_mitoBI <- RunPCA(object = seurC2_RNA_mitoBI, npcs = 50, verbose = FALSE)
seurC2_RNA_mitoBI <- RunUMAP(object = seurC2_RNA_mitoBI, reduction = "pca", dims = 1:50)
seurC2_RNA_mitoBI <- FindNeighbors(seurC2_RNA_mitoBI,reduction="pca",dims=1:50)
seurC2_RNA_mitoBI <- FindClusters(seurC2_RNA_mitoBI,resolution=0.8)

## pK Identification (no ground-truth)
sweep.res.list_seurC2_RNA_mitoBI <- paramSweep_v3(seurC2_RNA_mitoBI, PCs = 1:50, sct = FALSE)
sweep.stats_seurC2_RNA_mitoBI <- summarizeSweep(sweep.res.list_seurC2_RNA_mitoBI, GT = FALSE)
bcmvn_seurC2_RNA_mitoBI <- find.pK(sweep.stats_seurC2_RNA_mitoBI)
bcmvn_seurC2_RNA_mitoBI

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seurC2_RNA_mitoBI@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seurC2_RNA_mitoBI@meta.data$ClusteringResults
nExp_poi <- round(0.05*nrow(seurC2_RNA_mitoBI@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seurC2_RNA_mitoBI <- doubletFinder_v3(seurC2_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
#seurC2_RNA_mitoBI <- doubletFinder_v3(seurC2_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
DimPlot(seurC2_RNA_mitoBI, pt.size=0.1, group.by="DF.classifications_0.25_0.005_234")

## paste column to new place and delete old one
#seurC2_RNA_mitoBI$sample_doublets <- seurC2_RNA_mitoBI$DF.classifications_0.25_0.005_234
#seurC2_RNA_mitoBI$DF.classifications_0.25_0.005_234 <- NULL

##remove doublets and shrink object
seurC2_RNA_singlet<-subset(x=seurC2_RNA_mitoBI, subset = DF.classifications_0.25_0.005_234 == "Singlet")
seurC2_RNA_singlet_diet<-DietSeurat(seurC2_RNA_singlet)
seurC2_RNA_singlet_diet$DF.classifications_0.25_0.005_234 <-NULL
seurC2_RNA_singlet_diet$pANN_0.25_0.005_234<-NULL
seurC2_RNA_singlet_diet$seurat_clusters<-NULL
seurC2_RNA_singlet_diet$RNA_snn_res.0.8<-NULL
saveRDS(seurC2_RNA_singlet_diet,file="objects/seurC2_RNA.RDS")

DefaultAssay(object = seurC3_RNA_mitoBI) <- "RNA"
seurC3_RNA_mitoBI <- ScaleData(object = seurC3_RNA_mitoBI, verbose = FALSE)
seurC3_RNA_mitoBI <- FindVariableFeatures(seurC3_RNA_mitoBI)
seurC3_RNA_mitoBI <- RunPCA(object = seurC3_RNA_mitoBI, npcs = 50, verbose = FALSE)
seurC3_RNA_mitoBI <- RunUMAP(object = seurC3_RNA_mitoBI, reduction = "pca", dims = 1:50)
seurC3_RNA_mitoBI <- FindNeighbors(seurC3_RNA_mitoBI,reduction="pca",dims=1:50)
seurC3_RNA_mitoBI <- FindClusters(seurC3_RNA_mitoBI,resolution=0.8)

## pK Identification (no ground-truth)
sweep.res.list_seurC3_RNA_mitoBI <- paramSweep_v3(seurC3_RNA_mitoBI, PCs = 1:50, sct = FALSE)
sweep.stats_seurC3_RNA_mitoBI <- summarizeSweep(sweep.res.list_seurC3_RNA_mitoBI, GT = FALSE)
bcmvn_seurC3_RNA_mitoBI <- find.pK(sweep.stats_seurC3_RNA_mitoBI)
bcmvn_seurC3_RNA_mitoBI

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seurC3_RNA_mitoBI@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.05*nrow(seurC3_RNA_mitoBI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj


## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seurC3_RNA_mitoBI <- doubletFinder_v3(seurC3_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
#seurC3_RNA_mitoBI <- doubletFinder_v3(seurC3_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
DimPlot(seurC3_RNA_mitoBI, pt.size=0.1, group.by="DF.classifications_0.25_0.005_135")

## paste column to new place and delete old one
#seurC3_RNA_mitoBI$sample_doublets <- seurC3_RNA_mitoBI$DF.classifications_0.25_0.005_135
#seurC3_RNA_mitoBI$DF.classifications_0.25_0.005_135 <- NULL

#remove doublets and shrink object
seurC3_RNA_singlet<-subset(x=seurC3_RNA_mitoBI, subset = DF.classifications_0.25_0.005_135 == "Singlet")
seurC3_RNA_singlet_diet<-DietSeurat(seurC3_RNA_singlet)
seurC3_RNA_singlet_diet$DF.classifications_0.25_0.005_135 <-NULL
seurC3_RNA_singlet_diet$pANN_0.25_0.005_135<-NULL
seurC3_RNA_singlet_diet$seurat_clusters<-NULL
seurC3_RNA_singlet_diet$RNA_snn_res.0.8<-NULL
saveRDS(seurC3_RNA_singlet_diet,file="objects/seurC3_RNA.RDS")

DefaultAssay(object = seurC4_RNA_mitoBI) <- "RNA"
seurC4_RNA_mitoBI <- ScaleData(object = seurC4_RNA_mitoBI, verbose = FALSE)
seurC4_RNA_mitoBI <- FindVariableFeatures(seurC4_RNA_mitoBI)
seurC4_RNA_mitoBI <- RunPCA(object = seurC4_RNA_mitoBI, npcs = 50, verbose = FALSE)
seurC4_RNA_mitoBI <- RunUMAP(object = seurC4_RNA_mitoBI, reduction = "pca", dims = 1:50)
seurC4_RNA_mitoBI <- FindNeighbors(seurC4_RNA_mitoBI,reduction="pca",dims=1:50)
seurC4_RNA_mitoBI <- FindClusters(seurC4_RNA_mitoBI,resolution=0.8)

# pK Identification (no ground-truth)
sweep.res.list_seurC4_RNA_mitoBI <- paramSweep_v3(seurC4_RNA_mitoBI, PCs = 1:50, sct = FALSE)
sweep.stats_seurC4_RNA_mitoBI <- summarizeSweep(sweep.res.list_seurC4_RNA_mitoBI, GT = FALSE)
bcmvn_seurC4_RNA_mitoBI <- find.pK(sweep.stats_seurC4_RNA_mitoBI)
bcmvn_seurC4_RNA_mitoBI

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seurC4_RNA_mitoBI@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.05*nrow(seurC4_RNA_mitoBI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seurC4_RNA_mitoBI <- doubletFinder_v3(seurC4_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.04, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
#seurC4_RNA_mitoBI <- doC9_ALS_FTLD_RNA,file="Jan21_C9ALS_FTLD_RNA.RDS"ubletFinder_v3(seurC4_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
DimPlot(seurC4_RNA_mitoBI, pt.size=0.1, group.by="DF.classifications_0.25_0.05_257")

## paste column to new place and delete old one
#seurC4_RNA_mitoBI$sample_doublets <- seurC4_RNA_mitoBI$DF.classifications_0.25_0.04_257
#seurC4_RNA_mitoBI$DF.classifications_0.25_0.04_257 <- NULL

#remove doublets and shrink object
seurC4_RNA_singlet<-subset(x=seurC4_RNA_mitoBI, subset = DF.classifications_0.25_0.04_257 == "Singlet")
seurC4_RNA_singlet_diet<-DietSeurat(seurC4_RNA_singlet)
seurC4_RNA_singlet_diet$DF.classifications_0.25_0.04_257 <-NULL
seurC4_RNA_singlet_diet$pANN_0.25_0.04_257<-NULL
seurC4_RNA_singlet_diet$seurat_clusters<-NULL
seurC4_RNA_singlet_diet$RNA_snn_res.0.8<-NULL
saveRDS(seurC4_RNA_singlet_diet,file="objects/seurC4_RNA.RDS")

DefaultAssay(object = seurC6_RNA_mitoBI) <- "RNA"
seurC6_RNA_mitoBI <- ScaleData(object = seurC6_RNA_mitoBI, verbose = FALSE)
seurC6_RNA_mitoBI <- FindVariableFeatures(seurC6_RNA_mitoBI)
seurC6_RNA_mitoBI <- RunPCA(object = seurC6_RNA_mitoBI, npcs = 50, verbose = FALSE)
seurC6_RNA_mitoBI <- RunUMAP(object = seurC6_RNA_mitoBI, reduction = "pca", dims = 1:50)
seurC6_RNA_mitoBI <- FindNeighbors(seurC6_RNA_mitoBI,reduction="pca",dims=1:50)
seurC6_RNA_mitoBI <- FindClusters(seurC6_RNA_mitoBI,resolution=0.8)

# pK Identification (no ground-truth)
sweep.res.list_seurC6_RNA_mitoBI <- paramSweep_v3(seurC6_RNA_mitoBI, PCs = 1:50, sct = FALSE)
sweep.stats_seurC6_RNA_mitoBI <- summarizeSweep(sweep.res.list_seurC6_RNA_mitoBI, GT = FALSE)
bcmvn_seurC6_RNA_mitoBI <- find.pK(sweep.stats_seurC6_RNA_mitoBI)
bcmvn_seurC6_RNA_mitoBI

# Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seurC6_RNA_mitoBI@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.05*nrow(seurC6_RNA_mitoBI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seurC6_RNA_mitoBI <- doubletFinder_v3(seurC6_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.21, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
#seurC6_RNA_mitoBI <- doubletFinder_v3(seurC6_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
DimPlot(seurC6_RNA_mitoBI, pt.size=0.1, group.by="DF.classifications_0.25_0.21_433")

## paste column to new place and delete old one
#seurC6_RNA_mitoBI$sample_doublets <- seurC6_RNA_mitoBI$DF.classifications_0.25_0.21_433
#seurC6_RNA_mitoBI$DF.classifications_0.25_0.1_331 <- NULL

#remove doublets and shrink object
seurC6_RNA_singlet<-subset(x=seurC6_RNA_mitoBI, subset = DF.classifications_0.25_0.21_401 == "Singlet")
seurC6_RNA_singlet_diet<-DietSeurat(seurC6_RNA_singlet)
seurC6_RNA_singlet_diet$DF.classifications_0.25_0.21_401 <-NULL
seurC6_RNA_singlet_diet$pANN_0.25_0.21_401 <-NULL
seurC6_RNA_singlet_diet$seurat_clusters <-NULL
seurC6_RNA_singlet_diet$RNA_snn_res.0.8 <-NULL
saveRDS(seurC6_RNA_singlet_diet,file="objects/seurC6_RNA.RDS")