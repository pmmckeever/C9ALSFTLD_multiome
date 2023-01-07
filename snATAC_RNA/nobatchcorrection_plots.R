### sample DimPlots (correction vs. no correction)

library(Seurat)
library(Signac)
library(ggplot2)
library(patchwork)
library(scCustomize)
library(dittoSeq)
library(ggplotify)
library(ggpubr)
library(viridis)
set.seed(1234)

intRNA <- readRDS("objects/intRNA_final.RDS")
intATAC <- readRDS("objects/intATAC_final.RDS")

fragments<-CreateFragmentObject(path="fragments/fragments.tsv.gz",cells=colnames(intATAC))
Fragments(intATAC)<-NULL
Fragments(intATAC)<-fragments

#reorder samples for snRNA and snATAC
intRNA$sample<-factor(intRNA$sample, levels=c("CTRL1","CTRL2","CTRL3","CTRL4","CTRL5","CTRL6","CTRL7","C9noALSnoFTLD","C9ALSFTLD1","C9ALSFTLD2","C9ALSFTLD3","C9ALSFTLD4","C9ALSFTLD5","C9ALSFTLD6","C9ALSnoFTLD1","C9ALSnoFTLD2","C9ALSnoFTLD3","sALSnoFTLD1","sALSnoFTLD2","sALSnoFTLD3","sALSnoFTLD4","sALSnoFTLD5","sALSnoFTLD6","sALSnoFTLD7","sALSnoFTLD8"))
Idents(intRNA)<-"sample"
intATAC$sample<-factor(intATAC$sample, levels=c("CTRL1","CTRL2","CTRL3","CTRL4","C9noALSnoFTLD","C9ALSFTLD1","C9ALSFTLD2","C9ALSFTLD3","C9ALSFTLD5","C9ALSFTLD6","C9ALSFTLD7","C9ALSnoFTLD1","C9ALSnoFTLD2","C9ALSnoFTLD3","sALSnoFTLD1","sALSnoFTLD2","sALSnoFTLD3","sALSnoFTLD4","sALSnoFTLD5","sALSnoFTLD6","sALSnoFTLD7","sALSnoFTLD8"))
Idents(intATAC)<-"sample"

intRNA2 <- readRDS("objects/unint_RNA.RDS")

cells.use <- colnames(intRNA)
intRNA3 <- subset(intRNA2, cells = cells.use)
intRNA3
rm(intRNA2)
gc()
intRNA3<-FindVariableFeatures(intRNA3)
intRNA3<-ScaleData(intRNA3,vars.to.regress="percent.mt")
gc()
intRNA3<-RunPCA(intRNA3)
intRNA3<-RunUMAP(intRNA3, dims=1:30)
gc()
saveRDS(intRNA3,file="objects/RNA_nobatch.RDS")

intATAC2 <- readRDS("objects/unint_ATAC.RDS")
intATAC3 <- subset(intATAC2, cells=cells.use)
intATAC3 <- FindVariableFeatures(intATAC3)
intATAC3 <- RunTFIDF(intATAC3)
intATAC3 <- FindTopFeatures(intATAC3, min.cutoff = 50)
intATAC3 <- RunSVD(intATAC3)
intATAC3 <- RunUMAP(intATAC3, reduction = 'lsi', dims = 2:30)

saveRDS(intATAC3,file="objects/ATAC_nobatch.RDS")


sample_DimPlot_RNA<-DimPlot_scCustom(intRNA, group.by="sample",colors_use=pal_snRNA)+
  theme(legend.position="none")+ggtitle("snRNA by Sample")+labs(x="UMAP 1",y="UMAP 2")+
  mytheme+NoLegend()

nocorrection_sample_DimPlot_RNA<-DimPlot_scCustom(intRNA3, group.by="sample",colors_use=pal_snRNA)+
  theme(legend.position="right",legend.text=element_text(size=5), legend.key.size = unit(0.5, 'lines'))+
  guides(color = guide_legend(override.aes = list(size=2), ncol=2))+ggtitle("snRNA by Sample (no SCT integration)")+labs(x="UMAP 1",y="UMAP 2")+
  mytheme

gg_sample_DimPlot_RNA<-as.ggplot(sample_DimPlot_RNA)+theme(title=element_text(size=7))
gg_nocorrection_sample_DimPlot_RNA<-as.ggplot(nocorrection_sample_DimPlot_RNA)+theme(title=element_text(size=7))

sample_DimPlot_ATAC<-DimPlot_scCustom(intATAC, group.by="sample", colors_use=pal_snATAC)+
  theme(legend.position="none")+ggtitle("snATAC by Sample")+labs(x="UMAP 1",y="UMAP 2")+
  mytheme+NoLegend()

nocorrection_sample_DimPlot_ATAC<-DimPlot_scCustom(intATAC3, group.by="sample", colors_use=pal_snATAC, reduction="umap")+
  theme(legend.position="right",legend.text=element_text(size=5), legend.key.size = unit(0.5, 'lines'))+
  guides(color = guide_legend(override.aes = list(size=2), ncol=2))+ggtitle("snATAC by Sample (no Harmony)")+labs(x="UMAP 1",y="UMAP 2")+
  mytheme

gg_sample_DimPlot_ATAC<-as.ggplot(sample_DimPlot_ATAC)+theme(title=element_text(size=7))
gg_nocorrection_sample_DimPlot_ATAC<-as.ggplot(nocorrection_sample_DimPlot_ATAC)+theme(title=element_text(size=7))

