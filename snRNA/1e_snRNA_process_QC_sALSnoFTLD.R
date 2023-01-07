library(Seurat)
library(org.Hs.eg.db)
library(DoubletFinder)
set.seed(1234)

#sALSnoFTLD matrix location - COMPUTE CANADA CEDAR SERVER
input_from_10x1 <- "/ALS_snRNA/sALSnoFTLD1_snRNA/outs/filtered_feature_bc_matrix"
input_from_10x2 <- "/ALS_snRNA/sALSnoFTLD2_snRNA/outs/filtered_feature_bc_matrix"
input_from_10x3 <- "/ALS_snRNA/sALSnoFTLD3_snRNA/outs/filtered_feature_bc_matrix"
input_from_10x4 <- "/ALS_snRNA/sALSnoFTLD4_snRNA/outs/filtered_feature_bc_matrix"
input_from_10x5 <- "/ALS_snRNA/sALSnoFTLD5_snRNA/outs/filtered_feature_bc_matrix"
input_from_10x6 <- "/ALS_snRNA/sALSnoFTLD6_snRNA/outs/filtered_feature_bc_matrix"
input_from_10x7 <- "/ALS_snRNA/sALSnoFTLD7_snRNA/outs/filtered_feature_bc_matrix"
input_from_10x8 <- "/ALS_snRNA/sALSnoFTLD8_snRNA/outs/filtered_feature_bc_matrix"

#sALSnoFTLD cases read-in
seur1_RNA <- CreateSeuratObject(counts=Read10X(input_from_10x1), project="sALSnoFTLD1",
                          min.cells=1,min.features=1)
seur2_RNA <- CreateSeuratObject(counts=Read10X(input_from_10x2), project="sALSnoFTLD2",
                          min.cells=1,min.features=1)
seur3_RNA <- CreateSeuratObject(counts=Read10X(input_from_10x3), project="sALSnoFTLD3",
                            min.cells=1,min.features=1)
seur4_RNA <- CreateSeuratObject(counts=Read10X(input_from_10x4), project="sALSnoFTLD4",
                            min.cells=1,min.features=1)
seur5_RNA <- CreateSeuratObject(counts=Read10X(input_from_10x5), project="sALSnoFTLD5",
                          min.cells=1,min.features=1)
seur6_RNA <- CreateSeuratObject(counts=Read10X(input_from_10x6), project="sALSnoFTLD6",
                            min.cells=1,min.features=1)
seur7_RNA <- CreateSeuratObject(counts=Read10X(input_from_10x7), project="sALSnoFTLD7",
                          min.cells=1,min.features=1)
seur8_RNA <- CreateSeuratObject(counts=Read10X(input_from_10x8), project="sALSnoFTLD8",
                          min.cells=1,min.features=1)
show(seur1_RNA)
show(seur2_RNA)
show(seur3_RNA)
show(seur4_RNA)
show(seur5_RNA)
show(seur6_RNA)
show(seur7_RNA)
show(seur8_RNA)

seur1_RNA[["percent.mt"]] <- PercentageFeatureSet(seur1_RNA, pattern = "^MT-")
seur2_RNA[["percent.mt"]] <- PercentageFeatureSet(seur2_RNA, pattern = "^MT-")
seur3_RNA[["percent.mt"]] <- PercentageFeatureSet(seur3_RNA, pattern = "^MT-")
seur4_RNA[["percent.mt"]] <- PercentageFeatureSet(seur4_RNA, pattern = "^MT-")
seur5_RNA[["percent.mt"]] <- PercentageFeatureSet(seur5_RNA, pattern = "^MT-")
seur6_RNA[["percent.mt"]] <- PercentageFeatureSet(seur6_RNA, pattern = "^MT-")
seur7_RNA[["percent.mt"]] <- PercentageFeatureSet(seur7_RNA, pattern = "^MT-")
seur8_RNA[["percent.mt"]] <- PercentageFeatureSet(seur8_RNA, pattern = "^MT-")

seur1_RNA[["percent.rps"]] <- PercentageFeatureSet(seur1_RNA, pattern = "^RPS")
seur2_RNA[["percent.rps"]] <- PercentageFeatureSet(seur2_RNA, pattern = "^RPS")
seur3_RNA[["percent.rps"]] <- PercentageFeatureSet(seur3_RNA, pattern = "^RPS")
seur4_RNA[["percent.rps"]] <- PercentageFeatureSet(seur4_RNA, pattern = "^RPS")
seur5_RNA[["percent.rps"]] <- PercentageFeatureSet(seur5_RNA, pattern = "^RPS")
seur6_RNA[["percent.rps"]] <- PercentageFeatureSet(seur6_RNA, pattern = "^RPS")
seur7_RNA[["percent.rps"]] <- PercentageFeatureSet(seur7_RNA, pattern = "^RPS")
seur8_RNA[["percent.rps"]] <- PercentageFeatureSet(seur8_RNA, pattern = "^RPS")

seur1_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur1_RNA, pattern = "^RPL")
seur2_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur2_RNA, pattern = "^RPL")
seur3_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur3_RNA, pattern = "^RPL")
seur4_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur4_RNA, pattern = "^RPL")
seur5_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur5_RNA, pattern = "^RPL")
seur6_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur6_RNA, pattern = "^RPL")
seur7_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur7_RNA, pattern = "^RPL")
seur8_RNA[["percent.rpl"]] <- PercentageFeatureSet(seur8_RNA, pattern = "^RPL")


mito_gene_identifier <- "^MT-"
mads_thresh <- 3
hard_thresh <- 50
seur1_RNA <- PercentageFeatureSet(seur1_RNA, pattern = "^MT-", col.name = "percent.mt")
mito_thresh <- median(seur1_RNA$percent.mt) + mad(seur1_RNA$percent.mt) * mads_thresh
drop_mito <- seur1_RNA$percent.mt > mito_thresh | seur1_RNA$percent.mt > hard_thresh
seur1_RNA_mitoBI <- seur1_RNA[,!drop_mito]
show(seur1_RNA_mitoBI)
seur2_RNA <- PercentageFeatureSet(seur2_RNA, pattern = "^MT-", col.name = "percent.mt")
mito_thresh <- median(seur2_RNA$percent.mt) + mad(seur2_RNA$percent.mt) * mads_thresh
drop_mito <- seur2_RNA$percent.mt > mito_thresh | seur2_RNA$percent.mt > hard_thresh
seur2_RNA_mitoBI <- seur2_RNA[,!drop_mito]
show(seur2_RNA_mitoBI)
seur3_RNA <- PercentageFeatureSet(seur3_RNA, pattern = "^MT-", col.name = "percent.mt")
mito_thresh <- median(seur3_RNA$percent.mt) + mad(seur3_RNA$percent.mt) * mads_thresh
drop_mito <- seur3_RNA$percent.mt > mito_thresh | seur3_RNA$percent.mt > hard_thresh
seur3_RNA_mitoBI <- seur3_RNA[,!drop_mito]
show(seur3_RNA_mitoBI)
seur4_RNA <- PercentageFeatureSet(seur4_RNA, pattern = "^MT-", col.name = "percent.mt")
mito_thresh <- median(seur4_RNA$percent.mt) + mad(seur4_RNA$percent.mt) * mads_thresh
drop_mito <- seur4_RNA$percent.mt > mito_thresh | seur4_RNA$percent.mt > hard_thresh
seur4_RNA_mitoBI <- seur4_RNA[,!drop_mito]
show(seur4_RNA_mitoBI)
seur5_RNA <- PercentageFeatureSet(seur5_RNA, pattern = "^MT-", col.name = "percent.mt")
mito_thresh <- median(seur5_RNA$percent.mt) + mad(seur5_RNA$percent.mt) * mads_thresh
drop_mito <- seur5_RNA$percent.mt > mito_thresh | seur5_RNA$percent.mt > hard_thresh
seur5_RNA_mitoBI <- seur5_RNA[,!drop_mito]
show(seur5_RNA_mitoBI)
seur6_RNA <- PercentageFeatureSet(seur6_RNA, pattern = "^MT-", col.name = "percent.mt")
mito_thresh <- median(seur6_RNA$percent.mt) + mad(seur6_RNA$percent.mt) * mads_thresh
drop_mito <- seur6_RNA$percent.mt > mito_thresh | seur6_RNA$percent.mt > hard_thresh
seur6_RNA_mitoBI <- seur6_RNA[,!drop_mito]
show(seur6_RNA_mitoBI)
seur7_RNA <- PercentageFeatureSet(seur7_RNA, pattern = "^MT-", col.name = "percent.mt")
mito_thresh <- median(seur7_RNA$percent.mt) + mad(seur7_RNA$percent.mt) * mads_thresh
drop_mito <- seur7_RNA$percent.mt > mito_thresh | seur7_RNA$percent.mt > hard_thresh
seur7_RNA_mitoBI <- seur7_RNA[,!drop_mito]
show(seur7_RNA_mitoBI)
seur8_RNA <- PercentageFeatureSet(seur8_RNA, pattern = "^MT-", col.name = "percent.mt")
mito_thresh <- median(seur8_RNA$percent.mt) + mad(seur8_RNA$percent.mt) * mads_thresh
drop_mito <- seur8_RNA$percent.mt > mito_thresh | seur8_RNA$percent.mt > hard_thresh
seur8_RNA_mitoBI <- seur8_RNA[,!drop_mito]
show(seur8_RNA_mitoBI)

seur1_RNA_mitoBI <- subset(seur1_RNA_mitoBI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)
seur2_RNA_mitoBI <- subset(seur2_RNA_mitoBI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)
seur3_RNA_mitoBI <- subset(seur3_RNA_mitoBI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)
seur4_RNA_mitoBI <- subset(seur4_RNA_mitoBI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)
seur5_RNA_mitoBI <- subset(seur5_RNA_mitoBI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)
seur6_RNA_mitoBI <- subset(seur6_RNA_mitoBI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)
seur7_RNA_mitoBI <- subset(seur7_RNA_mitoBI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)
seur8_RNA_mitoBI <- subset(seur8_RNA_mitoBI, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)

show(seur1_RNA_mitoBI)
show(seur2_RNA_mitoBI)    
show(seur3_RNA_mitoBI)
show(seur4_RNA_mitoBI)
show(seur5_RNA_mitoBI)
show(seur6_RNA_mitoBI)
show(seur7_RNA_mitoBI)
show(seur8_RNA_mitoBI)

seur1_RNA_mitoBI$sample <- 'sALSnoFTLD1'
seur2_RNA_mitoBI$sample <- 'sALSnoFTLD2'
seur3_RNA_mitoBI$sample <- 'sALSnoFTLD3'
seur4_RNA_mitoBI$sample <- 'sALSnoFTLD4'
seur5_RNA_mitoBI$sample <- 'sALSnoFTLD5'
seur6_RNA_mitoBI$sample <- 'sALSnoFTLD6'
seur7_RNA_mitoBI$sample <- 'sALSnoFTLD7'
seur8_RNA_mitoBI$sample <- 'sALSnoFTLD8'


seur1_RNA_mitoBI<-NormalizeData(seur1_RNA_mitoBI)
seur2_RNA_mitoBI<-NormalizeData(seur2_RNA_mitoBI)    
seur3_RNA_mitoBI<-NormalizeData(seur3_RNA_mitoBI)
seur4_RNA_mitoBI<-NormalizeData(seur4_RNA_mitoBI)
seur5_RNA_mitoBI<-NormalizeData(seur5_RNA_mitoBI)
seur6_RNA_mitoBI<-NormalizeData(seur6_RNA_mitoBI)
seur7_RNA_mitoBI<-NormalizeData(seur7_RNA_mitoBI)
seur8_RNA_mitoBI<-NormalizeData(seur8_RNA_mitoBI)

seur1_RNA_mitoBI<-FindVariableFeatures(seur1_RNA_mitoBI)
seur2_RNA_mitoBI<-FindVariableFeatures(seur2_RNA_mitoBI)    
seur3_RNA_mitoBI<-FindVariableFeatures(seur3_RNA_mitoBI)
seur4_RNA_mitoBI<-FindVariableFeatures(seur4_RNA_mitoBI)
seur5_RNA_mitoBI<-FindVariableFeatures(seur5_RNA_mitoBI)
seur6_RNA_mitoBI<-FindVariableFeatures(seur6_RNA_mitoBI)
seur7_RNA_mitoBI<-FindVariableFeatures(seur7_RNA_mitoBI)
seur8_RNA_mitoBI<-FindVariableFeatures(seur8_RNA_mitoBI)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

seur1_RNA_mitoBI <- CellCycleScoring(seur1_RNA_mitoBI, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seur2_RNA_mitoBI <- CellCycleScoring(seur2_RNA_mitoBI, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seur3_RNA_mitoBI <- CellCycleScoring(seur3_RNA_mitoBI, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seur4_RNA_mitoBI <- CellCycleScoring(seur4_RNA_mitoBI, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seur5_RNA_mitoBI <- CellCycleScoring(seur5_RNA_mitoBI, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seur6_RNA_mitoBI <- CellCycleScoring(seur6_RNA_mitoBI, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seur7_RNA_mitoBI <- CellCycleScoring(seur7_RNA_mitoBI, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seur8_RNA_mitoBI <- CellCycleScoring(seur8_RNA_mitoBI, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DefaultAssay(object = seur1_RNA_mitoBI) <- "RNA"
seur1_RNA_mitoBI <- ScaleData(object = seur1_RNA_mitoBI, verbose = FALSE)
seur1_RNA_mitoBI <- FindVariableFeatures(seur1_RNA_mitoBI)
seur1_RNA_mitoBI <- RunPCA(object = seur1_RNA_mitoBI, npcs = 50, verbose = FALSE)
seur1_RNA_mitoBI <- RunUMAP(object = seur1_RNA_mitoBI, reduction = "pca", dims = 1:50)
seur1_RNA_mitoBI <- FindNeighbors(seur1_RNA_mitoBI,reduction="pca",dims=1:50)
seur1_RNA_mitoBI <- FindClusters(seur1_RNA_mitoBI,resolution=0.8)

## pK Identification (no ground-truth)
sweep.res.list_seur1_RNA_mitoBI <- paramSweep_v3(seur1_RNA_mitoBI, PCs = 1:50, sct = FALSE)
sweep.stats_seur1_RNA_mitoBI <- summarizeSweep(sweep.res.list_seur1_RNA_mitoBI, GT = FALSE)
bcmvn_seur1_RNA_mitoBI <- find.pK(sweep.stats_seur1_RNA_mitoBI)
bcmvn_seur1_RNA_mitoBI

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seur1_RNA_mitoBI@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.05*nrow(seur1_RNA_mitoBI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seur1_RNA_mitoBI <- doubletFinder_v3(seur1_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.02, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
#seur1_RNA_mitoBI <- doubletFinder_v3(seur1_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
DimPlot(seur1_RNA_mitoBI, pt.size=0.1, group.by="DF.classifications_0.25_0.02_97")

## paste column to new place and delete old one
seur1_RNA_mitoBI$sample_doublets <- seur1_RNA_mitoBI$DF.classifications_0.25_0.02_87
seur1_RNA_mitoBI$DF.classifications_0.25_0.03_87 <- NULL

##remove doublets and shrink object
seur1_RNA_singlet<-subset(x=seur1_RNA_mitoBI, subset = DF.classifications_0.25_0.02_87 == "Singlet")
seur1_RNA_singlet_diet<-DietSeurat(seur1_RNA_singlet)
seur1_RNA_singlet_diet$DF.classifications_0.25_0.02_87<-NULL
seur1_RNA_singlet_diet$pANN_0.25_0.02_87<-NULL
seur1_RNA_singlet_diet$seurat_clusters<-NULL
seur1_RNA_singlet_diet$RNA_snn_res.0.8<-NULL
saveRDS(seur1_RNA_singlet_diet,file="objects/seur1_RNA.RDS")

DefaultAssay(object = seur2_RNA_mitoBI) <- "RNA"
seur2_RNA_mitoBI <- ScaleData(object = seur2_RNA_mitoBI, verbose = FALSE)
seur2_RNA_mitoBI <- FindVariableFeatures(seur2_RNA_mitoBI)
seur2_RNA_mitoBI <- RunPCA(object = seur2_RNA_mitoBI, npcs = 50, verbose = FALSE)
seur2_RNA_mitoBI <- RunUMAP(object = seur2_RNA_mitoBI, reduction = "pca", dims = 1:50)
seur2_RNA_mitoBI <- FindNeighbors(seur2_RNA_mitoBI,reduction="pca",dims=1:50)
seur2_RNA_mitoBI <- FindClusters(seur2_RNA_mitoBI,resolution=0.8)

## pK Identification (no ground-truth)
sweep.res.list_seur2_RNA_mitoBI <- paramSweep_v3(seur2_RNA_mitoBI, PCs = 1:50, sct = FALSE)
sweep.stats_seur2_RNA_mitoBI <- summarizeSweep(sweep.res.list_seur2_RNA_mitoBI, GT = FALSE)
bcmvn_seur2_RNA_mitoBI <- find.pK(sweep.stats_seur2_RNA_mitoBI)
bcmvn_seur2_RNA_mitoBI

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seur2_RNA_mitoBI@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.05*nrow(seur2_RNA_mitoBI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seur2_RNA_mitoBI <- doubletFinder_v3(seur2_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
#seur2_RNA_mitoBI <- doubletFinder_v3(seur2_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
DimPlot(seur2_RNA_mitoBI, pt.size=0.1, group.by="DF.classifications_0.25_0.01_150")

## paste column to new place and delete old one
seur2_RNA_mitoBI$sample_doublets <- seur2_RNA_mitoBI$DF.classifications_0.01_150
seur2_RNA_mitoBI$DF.classifications_0.25_0.01_150 <- NULL

##remove doublets and shrink object
seur2_RNA_singlet<-subset(x=seur2_RNA_mitoBI, subset = DF.classifications_0.25_0.01_150 == "Singlet")
seur2_RNA_singlet_diet<-DietSeurat(seur2_RNA_singlet)
seur2_RNA_singlet_diet$DF.classifications_0.25_0.01_150<-NULL
seur2_RNA_singlet_diet$pANN_0.25_0.01_150<-NULL
seur2_RNA_singlet_diet$seurat_clusters<-NULL
seur2_RNA_singlet_diet$RNA_snn_res.0.8<-NULL
saveRDS(seur2_RNA_singlet_diet,file="objects/seur2_RNA.RDS")

DefaultAssay(object = seur3_RNA_mitoBI) <- "RNA"
seur3_RNA_mitoBI <- ScaleData(object = seur3_RNA_mitoBI, verbose = FALSE)
seur3_RNA_mitoBI <- FindVariableFeatures(seur3_RNA_mitoBI)
seur3_RNA_mitoBI <- RunPCA(object = seur3_RNA_mitoBI, npcs = 50, verbose = FALSE)
seur3_RNA_mitoBI <- RunUMAP(object = seur3_RNA_mitoBI, reduction = "pca", dims = 1:50)
seur3_RNA_mitoBI <- FindNeighbors(seur3_RNA_mitoBI,reduction="pca",dims=1:50)
seur3_RNA_mitoBI <- FindClusters(seur3_RNA_mitoBI,resolution=0.8)

## pK Identification (no ground-truth)
sweep.res.list_seur3_RNA_mitoBI <- paramSweep_v3(seur3_RNA_mitoBI, PCs = 1:50, sct = FALSE)
sweep.stats_seur3_RNA_mitoBI <- summarizeSweep(sweep.res.list_seur3_RNA_mitoBI, GT = FALSE)
bcmvn_seur3_RNA_mitoBI <- find.pK(sweep.stats_seur3_RNA_mitoBI)
bcmvn_seur3_RNA_mitoBI

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seur3_RNA_mitoBI@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.05*nrow(seur3_RNA_mitoBI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seur3_RNA_mitoBI <- doubletFinder_v3(seur3_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
#seur3_RNA_mitoBI <- doubletFinder_v3(seur3_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
DimPlot(seur3_RNA_mitoBI, pt.size=0.1, group.by="DF.classifications_0.25_0.005_252")

## paste column to new place and delete old one
seur3_RNA_mitoBI$sample_doublets <- seur3_RNA_mitoBI$DF.classifications_0.25_0.005_252
seur3_RNA_mitoBI$DF.classifications_0.25_0.005_252 <- NULL

##remove doublets and shrink object
seur3_RNA_singlet<-subset(x=seur3_RNA_mitoBI, subset = DF.classifications_0.25_0.005_252 == "Singlet")
seur3_RNA_singlet_diet<-DietSeurat(seur3_RNA_singlet)
seur3_RNA_singlet_diet$DF.classifications_0.25_0.005_252<-NULL
seur3_RNA_singlet_diet$pANN_0.25_0.005_252<-NULL
seur3_RNA_singlet_diet$seurat_clusters<-NULL
seur3_RNA_singlet_diet$RNA_snn_res.0.8<-NULL
saveRDS(seur3_RNA_singlet_diet,file="objects/seur3_RNA.RDS")

DefaultAssay(object = seur4_RNA_mitoBI) <- "RNA"
seur4_RNA_mitoBI <- ScaleData(object = seur4_RNA_mitoBI, verbose = FALSE)
seur4_RNA_mitoBI <- FindVariableFeatures(seur4_RNA_mitoBI)
seur4_RNA_mitoBI <- RunPCA(object = seur4_RNA_mitoBI, npcs = 50, verbose = FALSE)
seur4_RNA_mitoBI <- RunUMAP(object = seur4_RNA_mitoBI, reduction = "pca", dims = 1:50)
seur4_RNA_mitoBI <- FindNeighbors(seur4_RNA_mitoBI,reduction="pca",dims=1:50)
seur4_RNA_mitoBI <- FindClusters(seur4_RNA_mitoBI,resolution=0.8)

# pK Identification (no ground-truth)
sweep.res.list_seur4_RNA_mitoBI <- paramSweep_v3(seur4_RNA_mitoBI, PCs = 1:50, sct = FALSE)
sweep.stats_seur4_RNA_mitoBI <- summarizeSweep(sweep.res.list_seur4_RNA_mitoBI, GT = FALSE)
bcmvn_seur4_RNA_mitoBI <- find.pK(sweep.stats_seur4_RNA_mitoBI)
bcmvn_seur4_RNA_mitoBI

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seur4_RNA_mitoBI@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.05*nrow(seur4_RNA_mitoBI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seur4_RNA_mitoBI <- doubletFinder_v3(seur4_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.02, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
#seur4_RNA_mitoBI <- doubletFinder_v3(seur4_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
DimPlot(seur4_RNA_mitoBI, pt.size=0.1, group.by="DF.classifications_0.25_0.02_185")

## paste column to new place and delete old one
seur4_RNA_mitoBI$sample_doublets <- seur4_RNA_mitoBI$DF.classifications_0.25_0.02_185
seur4_RNA_mitoBI$DF.classifications_0.25_0.02_185 <- NULL

##remove doublets and shrink object
seur4_RNA_singlet<-subset(x=seur4_RNA_mitoBI, subset = DF.classifications_0.25_0.02_185 == "Singlet")
seur4_RNA_singlet_diet<-DietSeurat(seur4_RNA_singlet)
seur4_RNA_singlet_diet$DF.classifications_0.25_0.02_185<-NULL
seur4_RNA_singlet_diet$pANN_0.25_0.02_185<-NULL
seur4_RNA_singlet_diet$seurat_clusters<-NULL
seur4_RNA_singlet_diet$RNA_snn_res.0.8<-NULL
saveRDS(seur4_RNA_singlet_diet,file="objects/seur4_RNA.RDS")

DefaultAssay(object = seur5_RNA_mitoBI) <- "RNA"
seur5_RNA_mitoBI <- ScaleData(object = seur5_RNA_mitoBI, verbose = FALSE)
seur5_RNA_mitoBI <- FindVariableFeatures(seur5_RNA_mitoBI)
seur5_RNA_mitoBI <- RunPCA(object = seur5_RNA_mitoBI, npcs = 50, verbose = FALSE)
seur5_RNA_mitoBI <- RunUMAP(object = seur5_RNA_mitoBI, reduction = "pca", dims = 1:50)
seur5_RNA_mitoBI <- FindNeighbors(seur5_RNA_mitoBI,reduction="pca",dims=1:50)
seur5_RNA_mitoBI <- FindClusters(seur5_RNA_mitoBI,resolution=0.8)

# pK Identification (no ground-truth)
sweep.res.list_seur5_RNA_mitoBI <- paramSweep_v3(seur5_RNA_mitoBI, PCs = 1:50, sct = FALSE)
sweep.stats_seur5_RNA_mitoBI <- summarizeSweep(sweep.res.list_seur5_RNA_mitoBI, GT = FALSE)
bcmvn_seur5_RNA_mitoBI <- find.pK(sweep.stats_seur5_RNA_mitoBI)
bcmvn_seur5_RNA_mitoBI

# Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seur5_RNA_mitoBI@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.05*nrow(seur5_RNA_mitoBI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seur5_RNA_mitoBI <- doubletFinder_v3(seur5_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.19, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
#seur5_RNA_mitoBI <- doubletFinder_v3(seur5_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
DimPlot(seur5_RNA_mitoBI, pt.size=0.1, group.by="DF.classifications_0.25_0.19_263")

## paste column to new place and delete old one
seur5_RNA_mitoBI$sample_doublets <- seur5_RNA_mitoBI$DF.classifications_0.25_0.19_263
seur5_RNA_mitoBI$DF.classifications_0.25_0.19_263 <- NULL

##remove doublets and shrink object
seur5_RNA_singlet<-subset(x=seur5_RNA_mitoBI, subset = DF.classifications_0.25_0.19_263 == "Singlet")
seur5_RNA_singlet_diet<-DietSeurat(seur5_RNA_singlet)
seur5_RNA_singlet_diet$DF.classifications_0.25_0.19_263<-NULL
seur5_RNA_singlet_diet$pANN_0.25_0.19_263<-NULL
seur5_RNA_singlet_diet$seurat_clusters<-NULL
seur5_RNA_singlet_diet$RNA_snn_res.0.8<-NULL
saveRDS(seur5_RNA_singlet_diet,file="objects/seur5_RNA.RDS")

DefaultAssay(object = seur6_RNA_mitoBI) <- "RNA"
seur6_RNA_mitoBI <- ScaleData(object = seur6_RNA_mitoBI, verbose = FALSE)
seur6_RNA_mitoBI <- FindVariableFeatures(seur6_RNA_mitoBI)
seur6_RNA_mitoBI <- RunPCA(object = seur6_RNA_mitoBI, npcs = 50, verbose = FALSE)
seur6_RNA_mitoBI <- RunUMAP(object = seur6_RNA_mitoBI, reduction = "pca", dims = 1:50)
seur6_RNA_mitoBI <- FindNeighbors(seur6_RNA_mitoBI,reduction="pca",dims=1:50)
seur6_RNA_mitoBI <- FindClusters(seur6_RNA_mitoBI,resolution=0.8)

# pK Identification (no ground-truth)
sweep.res.list_seur6_RNA_mitoBI <- paramSweep_v3(seur6_RNA_mitoBI, PCs = 1:50, sct = FALSE)
sweep.stats_seur6_RNA_mitoBI <- summarizeSweep(sweep.res.list_seur6_RNA_mitoBI, GT = FALSE)
bcmvn_seur6_RNA_mitoBI <- find.pK(sweep.stats_seur6_RNA_mitoBI)
bcmvn_seur6_RNA_mitoBI

# Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seur6_RNA_mitoBI@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.05*nrow(seur6_RNA_mitoBI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj 

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seur6_RNA_mitoBI <- doubletFinder_v3(seur6_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.15, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
#seur6_RNA_mitoBI <- doubletFinder_v3(seur6_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
DimPlot(seur6_RNA_mitoBI, pt.size=0.1, group.by="DF.classifications_0.25_0.15_97")

## paste column to new place and delete old one
seur6_RNA_mitoBI$sample_doublets <- seur6_RNA_mitoBI$DF.classifications_0.25_0.15_97
seur6_RNA_mitoBI$DF.classifications_0.25_0.17_97 <- NULL

##remove doublets and shrink object
seur6_RNA_singlet<-subset(x=seur6_RNA_mitoBI, subset = DF.classifications_0.25_0.15_97 == "Singlet")
seur6_RNA_singlet_diet<-DietSeurat(seur6_RNA_singlet)
seur6_RNA_singlet_diet$DF.classifications_0.25_0.15_97<-NULL
seur6_RNA_singlet_diet$pANN_0.25_0.15_97<-NULL
seur6_RNA_singlet_diet$seurat_clusters<-NULL
seur6_RNA_singlet_diet$RNA_snn_res.0.8<-NULL
saveRDS(seur6_RNA_singlet_diet,file="objects/seur6_RNA.RDS")

DefaultAssay(object = seur7_RNA_mitoBI) <- "RNA"
seur7_RNA_mitoBI <- ScaleData(object = seur7_RNA_mitoBI, verbose = FALSE)
seur7_RNA_mitoBI <- FindVariableFeatures(seur7_RNA_mitoBI)
seur7_RNA_mitoBI <- RunPCA(object = seur7_RNA_mitoBI, npcs = 50, verbose = FALSE)
seur7_RNA_mitoBI <- RunUMAP(object = seur7_RNA_mitoBI, reduction = "pca", dims = 1:50)
seur7_RNA_mitoBI <- FindNeighbors(seur7_RNA_mitoBI,reduction="pca",dims=1:50)
seur7_RNA_mitoBI <- FindClusters(seur7_RNA_mitoBI,resolution=0.8)

## pK Identification (no ground-truth)
sweep.res.list_seur7_RNA_mitoBI <- paramSweep_v3(seur7_RNA_mitoBI, PCs = 1:50, sct = FALSE)
sweep.stats_seur7_RNA_mitoBI <- summarizeSweep(sweep.res.list_seur7_RNA_mitoBI, GT = FALSE)
bcmvn_seur7_RNA_mitoBI <- find.pK(sweep.stats_seur7_RNA_mitoBI)
bcmvn_seur7_RNA_mitoBI

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seur7_RNA_mitoBI@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)         
nExp_poi <- round(0.05*nrow(seur7_RNA_mitoBI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seur7_RNA_mitoBI <- doubletFinder_v3(seur7_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
#seur7_RNA_mitoBI <- doubletFinder_v3(seur7_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
DimPlot(seur7_RNA_mitoBI, pt.size=0.1, group.by="DF.classifications_0.25_0.005_150")

## paste column to new place and delete old one
seur7_RNA_mitoBI$sample_doublets <- seur7_RNA_mitoBI$DF.classifications_0.25_0.005_150
seur7_RNA_mitoBI$DF.classifications_0.25_0.005_150 <- NULL

#remove doublets and shrink object
seur7_RNA_singlet<-subset(x=seur7_RNA_mitoBI, subset = DF.classifications_0.25_0.005_150 == "Singlet")
seur7_RNA_singlet_diet<-DietSeurat(seur7_RNA_singlet)
seur7_RNA_singlet_diet$DF.classifications_0.25_0.005_150 <-NULL
seur7_RNA_singlet_diet$pANN_0.25_0.005_150<-NULL
seur7_RNA_singlet_diet$seurat_clusters<-NULL
seur7_RNA_singlet_diet$RNA_snn_res.0.8<-NULL
saveRDS(seur7_RNA_singlet_diet,file="objects/seur7_RNA.RDS")

DefaultAssay(object = seur8_RNA_mitoBI) <- "RNA"
seur8_RNA_mitoBI <- ScaleData(object = seur8_RNA_mitoBI, verbose = FALSE)
seur8_RNA_mitoBI <- FindVariableFeatures(seur8_RNA_mitoBI)
seur8_RNA_mitoBI <- RunPCA(object = seur8_RNA_mitoBI, npcs = 50, verbose = FALSE)
seur8_RNA_mitoBI <- RunUMAP(object = seur8_RNA_mitoBI, reduction = "pca", dims = 1:50)
seur8_RNA_mitoBI <- FindNeighbors(seur8_RNA_mitoBI,reduction="pca",dims=1:50)
seur8_RNA_mitoBI <- FindClusters(seur8_RNA_mitoBI,resolution=0.8)

## pK Identification (no ground-truth)
sweep.res.list_seur8_RNA_mitoBI <- paramSweep_v3(seur8_RNA_mitoBI, PCs = 1:50, sct = FALSE)
sweep.stats_seur8_RNA_mitoBI <- summarizeSweep(sweep.res.list_seur8_RNA_mitoBI, GT = FALSE)
bcmvn_seur8_RNA_mitoBI <- find.pK(sweep.stats_seur8_RNA_mitoBI)
bcmvn_seur8_RNA_mitoBI

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- seur8_RNA_mitoBI@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)         
nExp_poi <- round(0.05*nrow(seur8_RNA_mitoBI@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
nExp_poi.adj

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seur8_RNA_mitoBI <- doubletFinder_v3(seur8_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.11, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
#seur8_RNA_mitoBI <- doubletFinder_v3(seur8_RNA_mitoBI, PCs = 1:50, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
DimPlot(seur8_RNA_mitoBI, pt.size=0.1, group.by="DF.classifications_0.25_0.11_211")

## paste column to new place and delete old one
seur8_RNA_mitoBI$sample_doublets <- seur8_RNA_mitoBI$DF.classifications_0.25_0.11_211
seur8_RNA_mitoBI$DF.classifications_0.25_0.11_211 <- NULL

#remove doublets and shrink object
seur8_RNA_singlet<-subset(x=seur8_RNA_mitoBI, subset = DF.classifications_0.25_0.11_211 == "Singlet")
seur8_RNA_singlet_diet<-DietSeurat(seur8_RNA_singlet)
seur8_RNA_singlet_diet$DF.classifications_0.25_0.11_211 <-NULL
seur8_RNA_singlet_diet$pANN_0.25_0.11_211<-NULL
seur8_RNA_singlet_diet$seurat_clusters<-NULL
seur8_RNA_singlet_diet$RNA_snn_res.0.8<-NULL
saveRDS(seur8_RNA_singlet_diet,file="objects/seur8_RNA.RDS")