library(Seurat)
library(Signac)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(BiocParallel)
register(MulticoreParam(6)) 
set.seed(1234)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- 'hg38'
pfm <- getMatrixSet(x = JASPAR2020,opts = list(species = 9606, all_versions = FALSE))
pfm
ATAC_links <- readRDS("objects/ATAC_links.RDS")

intMICRO <- readRDS("objects/ATAC_microglia.RDS")
DefaultAssay(intMICRO) <- "peaks"
Annotation(intMICRO) <- annotations
DefaultAssay(intMICRO)<-'chromvar'
rownames(intMICRO)
motif_IDs<-rownames(intMICRO)
head(motif_IDs)
DefaultAssay(intMICRO)<-"peaks"
motif_names<-ConvertMotifID(intMICRO,id=motif_IDs)
names(pfm)<-motif_names
pfm

### MICROGLIA
intMICRO2<-DietSeurat(intMICRO,counts=T,data=T,scale.data=T,assays=c("peaks","RNA","chromvar"),dimreducs=c("lsi", "umap", "harmony"))

#Seurat not transfering harmony and lsi, so move embeddings manually
intMICRO2@reductions$harmony <- intMICRO@reductions$harmony
intMICRO2@reductions$lsi <- intMICRO@reductions$lsi

intMICRO2 <- AddMotifs(intMICRO2,genome=BSgenome.Hsapiens.UCSC.hg38, pfm=pfm)
gc()

intMICRO2 <- RunChromVAR(object = intMICRO2,genome = BSgenome.Hsapiens.UCSC.hg38)
intMICRO2

intMICRO2<-RegionStats(intMICRO2,genome=BSgenome.Hsapiens.UCSC.hg38)
intMICRO2

Links(intMICRO2)<-ATAC_links

saveRDS(intMICRO2,file="objects/ATAC_microglia.RDS")

DefaultAssay(intMICRO2)<-"chromvar"

intMICRO2$diagnosis<-intMICRO2$group
Idents(intMICRO2)<-"diagnosis"

avg_chromvar_micro<-AverageExpression(intMICRO2,assays="chromvar",return.seurat=T)
avg_chromvar_micro$diagnosis<-Idents(avg_chromvar_micro)
saveRDS(avg_chromvar_micro, file="objects/avg_chromvar_micro.RDS")

### OPCs
intOPC <- readRDS("objects/ATAC_OPC.RDS")

DefaultAssay(intOPC) <- "peaks"
Annotation(intOPC) <- annotations

intOPC2<-DietSeurat(intOPC,counts=T,data=T,scale.data=T,assays=c("peaks","RNA","chromvar"),dimreducs=c("lsi", "umap", "harmony"))

#Seurat not transfering harmony and lsi, so move embeddings manually
intOPC2@reductions$harmony <- intOPC@reductions$harmony
intOPC2@reductions$lsi <- intOPC@reductions$lsi

intOPC2 <- AddMotifs(intOPC2,genome=BSgenome.Hsapiens.UCSC.hg38, pfm=pfm)
gc()

intOPC2 <- RunChromVAR(object = intOPC2,genome = BSgenome.Hsapiens.UCSC.hg38)
intOPC2

intOPC2<-RegionStats(intOPC2,genome=BSgenome.Hsapiens.UCSC.hg38)
intOPC2

Links(intOPC2)<-ATAC_links

saveRDS(intOPC2,file="objects/ATAC_OPC.RDS")

DefaultAssay(intOPC2)<-"chromvar"

intOPC2$diagnosis<-intOPC2$group
Idents(intOPC2)<-"diagnosis"

avg_chromvar_OPC<-AverageExpression(intOPC2,assays="chromvar",return.seurat=T)
avg_chromvar_OPC$diagnosis<-Idents(avg_chromvar_OPC)
saveRDS(avg_chromvar_OPC,file="objects/avg_chromvar_OPC.RDS")

### ASTROCYTES

intAST <- readRDS("objects/ATAC_astrocytes.RDS")

DefaultAssay(intAST) <- "peaks"
Annotation(intAST) <- annotations

intAST2<-DietSeurat(intAST,counts=T,data=T,scale.data=T,assays=c("peaks","RNA","chromvar"),dimreducs=c("lsi", "umap", "harmony"))

#Seurat not transfering harmony and lsi, so move embeddings manually
intAST2@reductions$harmony <- intAST@reductions$harmony
intAST2@reductions$lsi <- intAST@reductions$lsi

intAST2 <- AddMotifs(intAST2,genome=BSgenome.Hsapiens.UCSC.hg38, pfm=pfm)
gc()

intAST2 <- RunChromVAR(object = intAST2,genome = BSgenome.Hsapiens.UCSC.hg38)
intAST2

intAST2<-RegionStats(intAST2,genome=BSgenome.Hsapiens.UCSC.hg38)
intAST2

Links(intAST2)<-ATAC_links

saveRDS(intAST2,file="objects/ATAC_astrocytes.RDS")

DefaultAssay(intAST2)<-"chromvar"

intAST2$diagnosis<-intAST$group
Idents(intAST2)<-"diagnosis"

avg_chromvar_astro<-AverageExpression(intAST2,assays="chromvar",return.seurat=T)
avg_chromvar_astro$diagnosis<-Idents(avg_chromvar_astro)
saveRDS(avg_chromvar_astro, file="objects/avg_chromvar_astro.RDS")


### EXCITATORY
intEXC <- readRDS("objects/ATAC_excitatory.RDS")
### OLIGODENDROCYTES

intOLIGO <- readRDS("objects/ATAC_oligo.RDS")

DefaultAssay(intOLIGO) <- "peaks"
Annotation(intOLIGO) <- annotations

intOLIGO2<-DietSeurat(intOLIGO,counts=T,data=T,scale.data=T,assays=c("peaks","RNA","chromvar"),dimreducs=c("lsi", "umap", "harmony"))

#Seurat not transfering harmony and lsi, so move embeddings manually
intOLIGO2@reductions$harmony <- intOLIGO@reductions$harmony
intOLIGO2@reductions$lsi <- intOLIGO@reductions$lsi

intOLIGO2 <- AddMotifs(intOLIGO2,genome=BSgenome.Hsapiens.UCSC.hg38, pfm=pfm)
gc()

intOLIGO2 <- RunChromVAR(object = intOLIGO2,genome = BSgenome.Hsapiens.UCSC.hg38)
intOLIGO2

intOLIGO2<-RegionStats(intOLIGO2,genome=BSgenome.Hsapiens.UCSC.hg38)
intOLIGO2

Links(intOLIGO2)<-ATAC_links

saveRDS(intOLIGO2,file="objects/ATAC_oligo.RDS")

DefaultAssay(intOLIGO2)<-"chromvar"

intOLIGO2$diagnosis<-intOLIGO2$group
Idents(intOLIGO2)<-"diagnosis"

avg_chromvar_oligo<-AverageExpression(intOLIGO2,assays="chromvar",return.seurat=T)
avg_chromvar_oligo$diagnosis<-Idents(avg_chromvar_oligo)
saveRDS(avg_chromvar_oligo,file="objects/avg_chromvar_oligo.RDS")


DefaultAssay(intEXC) <- "peaks"
Annotation(intEXC) <- annotations

intEXC2<-DietSeurat(intEXC,counts=T,data=T,scale.data=T,assays=c("peaks","RNA","chromvar"),dimreducs=c("lsi", "umap", "harmony"))

#Seurat not transfering harmony and lsi, so move embeddings manually
intEXC2@reductions$harmony <- intEXC@reductions$harmony
intEXC2@reductions$lsi <- intEXC@reductions$lsi

intEXC2 <- AddMotifs(intEXC2,genome=BSgenome.Hsapiens.UCSC.hg38, pfm=pfm)
gc()

intEXC2 <- RunChromVAR(object = intEXC2,genome = BSgenome.Hsapiens.UCSC.hg38)
intEXC2

intEXC2<-RegionStats(intEXC2,genome=BSgenome.Hsapiens.UCSC.hg38)
intEXC2

Links(intEXC2)<-ATAC_links

saveRDS(intEXC2,file="objects/ATAC_excitatory.RDS")

DefaultAssay(intEXC2)<-"chromvar"

intEXC2$diagnosis<-intEXC$group
Idents(intEXC2)<-"diagnosis"

avg_chromvar_ex<-AverageExpression(intEXC2,assays="chromvar",return.seurat=T)
avg_chromvar_ex$diagnosis<-Idents(avg_chromvar_ex)
saveRDS(avg_chromvar_ex, file="objects/avg_chromvar_excitatory.RDS")

### INHIBITORY

intINH <- readRDS("objects/ATAC_inhibitory.RDS")

DefaultAssay(intINH) <- "peaks"
Annotation(intINH) <- annotations

intINH2<-DietSeurat(intINH,counts=T,data=T,scale.data=T,assays=c("peaks","RNA","chromvar"),dimreducs=c("lsi", "umap", "harmony"))

#Seurat not transfering harmony and lsi, so move embeddings manually
intINH2@reductions$harmony <- intINH@reductions$harmony
intINH2@reductions$lsi <- intINH@reductions$lsi

intINH2 <- AddMotifs(intINH2,genome=BSgenome.Hsapiens.UCSC.hg38, pfm=pfm)
gc()

intINH2 <- RunChromVAR(object = intINH2,genome = BSgenome.Hsapiens.UCSC.hg38)
intINH2

intINH2<-RegionStats(intINH2,genome=BSgenome.Hsapiens.UCSC.hg38)
intINH2

Links(intINH2)<-ATAC_links

saveRDS(intINH2,file="objects/ATAC_inhibitory.RDS")

DefaultAssay(intINH2)<-"chromvar"

intINH2$diagnosis<-intINH$group
Idents(intINH2)<-"diagnosis"

avg_chromvar_in<-AverageExpression(intINH2,assays="chromvar",return.seurat=T)
avg_chromvar_in$diagnosis<-Idents(avg_chromvar_in)
saveRDS(avg_chromvar_in, file="objects/avg_chromvar_inhibitory.RDS")

### OLIGODENDROCYTES

intOLIGO <- readRDS("objects/ATAC_oligo.RDS")

DefaultAssay(intOLIGO) <- "peaks"
Annotation(intOLIGO) <- annotations

intOLIGO2<-DietSeurat(intOLIGO,counts=T,data=T,scale.data=T,assays=c("peaks","RNA","chromvar"),dimreducs=c("lsi", "umap", "harmony"))

#Seurat not transfering harmony and lsi, so move embeddings manually
intOLIGO2@reductions$harmony <- intOLIGO@reductions$harmony
intOLIGO2@reductions$lsi <- intOLIGO@reductions$lsi

intOLIGO2 <- AddMotifs(intOLIGO2,genome=BSgenome.Hsapiens.UCSC.hg38, pfm=pfm)
gc()

intOLIGO2 <- RunChromVAR(object = intOLIGO2,genome = BSgenome.Hsapiens.UCSC.hg38)
intOLIGO2

intOLIGO2<-RegionStats(intOLIGO2,genome=BSgenome.Hsapiens.UCSC.hg38)
intOLIGO2

Links(intOLIGO2)<-ATAC_links

saveRDS(intOLIGO2,file="objects/ATAC_oligo.RDS")

DefaultAssay(intOLIGO2)<-"chromvar"

intOLIGO2$diagnosis<-intOLIGO2$group
Idents(intOLIGO2)<-"diagnosis"

avg_chromvar_oligo<-AverageExpression(intOLIGO2,assays="chromvar",return.seurat=T)
avg_chromvar_oligo$diagnosis<-Idents(avg_chromvar_oligo)
saveRDS(avg_chromvar_oligo,file="objects/avg_chromvar_oligo.RDS")

#################### WHOLE DATASET

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- 'hg38'
pfm <- getMatrixSet(x = JASPAR2020,opts = list(species = 9606, all_versions = FALSE))
pfm
ATAC_links <- readRDS("objects/ATAC_links.RDS")

intATAC <- readRDS("objects/intATAC_final.RDS")
DefaultAssay(intATAC) <- "peaks"
Annotation(intATAC) <- annotations
DefaultAssay(intATAC)<-'chromvar'
rownames(intATAC)
motif_IDs<-rownames(intATAC)
head(motif_IDs)
DefaultAssay(intATAC)<-"peaks"
motif_names<-ConvertMotifID(intATAC,id=motif_IDs)
names(pfm)<-motif_names
pfm
Annotation(intATAC) <- annotations

intATAC2<-DietSeurat(intATAC,counts=T,data=T,scale.data=T,assays=c("peaks","RNA","chromvar"),dimreducs=c("lsi", "umap", "harmony"))

#Seurat not transfering harmony and lsi, so move embeddings manually
intATAC2@reductions$harmony <- intATAC@reductions$harmony
intATAC2@reductions$lsi <- intATAC@reductions$lsi

intATAC2 <- AddMotifs(intATAC2,genome=BSgenome.Hsapiens.UCSC.hg38, pfm=pfm)
gc()

intATAC2 <- RunChromVAR(object = intATAC2,genome = BSgenome.Hsapiens.UCSC.hg38)
intATAC2

intATAC2<-RegionStats(intATAC2,genome=BSgenome.Hsapiens.UCSC.hg38)
intATAC2

Links(intATAC2)<-ATAC_links

saveRDS(intATAC2,file="objects/intATAC_final.RDS")

DefaultAssay(intATAC2)<-"chromvar"

intATAC2$diagnosis<-intATAC2$group
Idents(intATAC2)<-"diagnosis"

avg_chromvar_atac<-AverageExpression(intATAC2,assays="chromvar",return.seurat=T)
avg_chromvar_atac$diagnosis<-Idents(avg_chromvar_atac)
saveRDS(avg_chromvar_atac,file="objects/avg_chromvar_ATAC.RDS")

RENAME MOTIF NAMES IN FINDMARKERS FILE
motif_names<-ConvertMotifID(Motifs(intATAC),id=celltype_motif_markers$gene,assay='peaks')
celltype_motif_markers$gene<-motif_names

colnames<-rownames()

motif_names<-ConvertMotifID(Motifs(intATAC),id=)
names(pfm)<-motif_names
