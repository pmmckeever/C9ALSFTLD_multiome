
library(Seurat)
set.seed(1234)
intRNA <- readRDS("objects/intRNA_V3only.RDS")

celltypes <-c ("Oligodendrocytes","OPC","AST-PP","AST-FB","Microglia","L23","L4","L56","L56-CC","IN-PV","IN-SST","IN-SV2C","IN-VIP")
Idents(intRNA)<-'cellsubtype'
combined_df <- data.frame()
for(i in 1:length(celltypes)){
  
  cur_celltype <- celltypes[[i]]
  print(cur_celltype)
  
  cur_seurat_obj <- subset(intRNA, cellsubtype == cur_celltype)
  Idents(cur_seurat_obj) <- cur_seurat_obj$diagnosis
  print(unique(Idents(cur_seurat_obj)))
  
  markers <- FindMarkers(object = cur_seurat_obj, 
                         ident.1=c("C9ALSFTLD"),
                         ident.2=c("control"),
                         test.use = "LR",
                         latent.vars = c("sex","nCount_RNA"),
                         logfc.threshold = 0.1,
                         min.pct = 0.3,
                         min.diff.pct = -Inf,
                         verbose = TRUE,
                         random.seed = 1,
                         pseudocount.use = 1,
                         assay='RNA')
  markers$diff <- markers$pct.1 - markers$pct.2
  markers$gene <- rownames(markers)
  markers$cellsubtype <- cur_celltype
  
  combined_df <- rbind(combined_df, markers)
}
cellsubtype.C9ALSFTLD.markers <- combined_df

save(cellsubtype.C9ALSFTLD.markers,file="cellsubtype_C9ALSFTLD_markers_LR_RNA_V3only_wRP.rda")
write.csv(cellsubtype.C9ALSFTLD.markers, file='cellsubtype_C9ALSFTLD_markers_LR_RNA_V3only_wRP.csv', quote=F)

rm(combined_df)

Idents(intRNA)<-'cellsubtype'
combined_df <- data.frame()
for(i in 1:length(celltypes)){
  
  cur_celltype <- celltypes[[i]]
  print(cur_celltype)
  
  cur_seurat_obj <- subset(intRNA, cellsubtype == cur_celltype)
  Idents(cur_seurat_obj) <- cur_seurat_obj$diagnosis
  print(unique(Idents(cur_seurat_obj)))
  
  markers <- FindMarkers(object = cur_seurat_obj, 
                         ident.1=c("C9ALSnoFTLD"),
                         ident.2=c("control"),
                         test.use = "LR",
                         latent.vars = c("sex","nCount_RNA"),
                         logfc.threshold = 0.1,
                         min.pct = 0.3,
                         min.diff.pct = -Inf,
                         verbose = TRUE,
                         random.seed = 1,
                         pseudocount.use = 1,
                         assay='RNA')
  markers$diff <- markers$pct.1 - markers$pct.2
  markers$gene <- rownames(markers)
  markers$cellsubtype <- cur_celltype
  
  combined_df <- rbind(combined_df, markers)
}
cellsubtype.C9ALSnoFTLD.markers <- combined_df

save(cellsubtype.C9ALSnoFTLD.markers,file="cellsubtype_C9ALSnoFTLD_markers_LR_RNA_V3only_wRP.rda")
write.csv(cellsubtype.C9ALSnoFTLD.markers, file='cellsubtype_C9ALSnoFTLD_markers_LR_RNA_V3only_wRP.csv', quote=F)

rm(combined_df)

Idents(intRNA)<-'cellsubtype'
combined_df <- data.frame()
for(i in 1:length(celltypes)){
  
  cur_celltype <- celltypes[[i]]
  print(cur_celltype)
  
  cur_seurat_obj <- subset(intRNA, cellsubtype == cur_celltype)
  Idents(cur_seurat_obj) <- cur_seurat_obj$diagnosis
  print(unique(Idents(cur_seurat_obj)))
  
  markers <- FindMarkers(object = cur_seurat_obj, 
                         ident.1=c("sALSnoFTLD"),
                         ident.2=c("control"),
                         test.use = "LR",
                         latent.vars = c("sex","nCount_RNA"),
                         logfc.threshold = 0.1,
                         min.pct = 0.3,
                         min.diff.pct = -Inf,
                         verbose = TRUE,
                         random.seed = 1,
                         pseudocount.use = 1,
                         assay='RNA')
  markers$diff <- markers$pct.1 - markers$pct.2
  markers$gene <- rownames(markers)
  markers$cellsubtype <- cur_celltype
  
  combined_df <- rbind(combined_df, markers)
}
cellsubtype.sALSnoFTLD.markers <- combined_df

save(cellsubtype.sALSnoFTLD.markers,file="cellsubtype_sALSnoFTLD_markers_LR_RNA_V3only_wRP.rda")
write.csv(cellsubtype.sALSnoFTLD.markers, file='cellsubtype_sALSnoFTLD_markers_LR_RNA_V3only_wRP.csv', quote=F)

rm(combined_df)

Idents(intRNA)<-'cellsubtype'
combined_df <- data.frame()
for(i in 1:length(celltypes)){
  
  cur_celltype <- celltypes[[i]]
  print(cur_celltype)
  
  cur_seurat_obj <- subset(intRNA, cellsubtype == cur_celltype)
  Idents(cur_seurat_obj) <- cur_seurat_obj$diagnosis
  print(unique(Idents(cur_seurat_obj)))
  
  markers <- FindMarkers(object = cur_seurat_obj, 
                         ident.1=c("C9ALSFTLD","C9ALSnoFTLD"),
                         ident.2=c("control"),
                         test.use = "LR",
                         latent.vars = c("sex","nCount_RNA"),
                         logfc.threshold = 0.1,
                         min.pct = 0.3,
                         min.diff.pct = -Inf,
                         verbose = TRUE,
                         random.seed = 1,
                         pseudocount.use = 1,
                         assay='RNA')
  markers$diff <- markers$pct.1 - markers$pct.2
  markers$gene <- rownames(markers)
  markers$cellsubtype <- cur_celltype
  
  combined_df <- rbind(combined_df, markers)
}
cellsubtype.C9together.markers <- combined_df

save(cellsubtype.C9together.markers,file="cellsubtype_C9together_markers_LR_RNA_V3only_wRP.rda")
write.csv(cellsubtype.C9together.markers, file='cellsubtype_C9together_markers_LR_RNA_V3only_wRP.csv', quote=F)

##########################################
##########################################
##########################################

rm(combined_df)

celltypes <- c("Oligodendrocytes","OPC","Astrocytes","Microglia","Excitatory","Inhibitory")
Idents(intRNA)<-'celltype'
combined_df <- data.frame()
for(i in 1:length(celltypes)){
  
  cur_celltype <- celltypes[[i]]
  print(cur_celltype)
  
  cur_seurat_obj <- subset(intRNA, celltype == cur_celltype)
  Idents(cur_seurat_obj) <- cur_seurat_obj$diagnosis
  print(unique(Idents(cur_seurat_obj)))
  
  markers <- FindMarkers(object = cur_seurat_obj, 
                         ident.1=c('C9ALSFTLD'),
                         ident.2=c('control'),
                         test.use = "LR",
                         latent.vars = c("sex","nCount_RNA"),
                         logfc.threshold = 0.1,
                         min.pct = 0.3,
                         min.diff.pct = -Inf,
                         verbose = TRUE,
                         random.seed = 1,
                         pseudocount.use = 1,
                         assay='RNA')
  markers$diff <- markers$pct.1 - markers$pct.2
  markers$gene <- rownames(markers)
  markers$celltype <- cur_celltype
  
  combined_df <- rbind(combined_df, markers)
}
celltype.C9ALSFTLD.markers <- combined_df

save(celltype.C9ALSFTLD.markers,file="celltype_C9ALSFTLD_markers_LR_RNA.rda")
write.csv(celltype.C9ALSFTLD.markers, file='celltype_C9ALSFTLD_markers_LR_RNA.csv', quote=F)

combined_df2 <- data.frame()
for(i in 1:length(celltypes)){
  
  cur_celltype <- celltypes[[i]]
  print(cur_celltype)
  
  cur_seurat_obj <- subset(intRNA, celltype == cur_celltype)
  Idents(cur_seurat_obj) <- cur_seurat_obj$diagnosis
  print(unique(Idents(cur_seurat_obj)))
  
  markers <- FindMarkers(object = cur_seurat_obj, 
                         ident.1=c('C9ALSnoFTLD'),
                         ident.2=c('control'),
                         test.use = "LR",
                         latent.vars = c("sex","nCount_RNA"),
                         logfc.threshold = 0.1,
                         min.pct = 0.3,
                         min.diff.pct = -Inf,
                         verbose = TRUE,
                         random.seed = 1,
                         pseudocount.use = 1,
                         assay='RNA')
  markers$diff <- markers$pct.1 - markers$pct.2
  markers$gene <- rownames(markers)
  markers$celltype <- cur_celltype
  
  combined_df2 <- rbind(combined_df2, markers)
}
celltype.C9ALSnoFTLD.markers <- combined_df2

save(celltype.C9ALSnoFTLD.markers,file="celltype_C9ALSnoFTLD_markers_LR_RNA.rda")
write.csv(celltype.C9ALSnoFTLD.markers, file='celltype_C9ALSnoFTLD_markers_LR_RNA.csv', quote=F)

combined_df3 <- data.frame()
for(i in 1:length(celltypes)){
  
  cur_celltype <- celltypes[[i]]
  print(cur_celltype)
  
  cur_seurat_obj <- subset(intRNA, celltype == cur_celltype)
  Idents(cur_seurat_obj) <- cur_seurat_obj$diagnosis
  print(unique(Idents(cur_seurat_obj)))
  
  markers <- FindMarkers(object = cur_seurat_obj, 
                         ident.1=c('sALSnoFTLD'),
                         ident.2=c('control'),
                         test.use = "LR",
                         latent.vars = c("sex","nCount_RNA"),
                         logfc.threshold = 0.1,
                         min.pct = 0.3,
                         min.diff.pct = -Inf,
                         verbose = TRUE,
                         random.seed = 1,
                         pseudocount.use = 1,
                         assay='RNA')
  markers$diff <- markers$pct.1 - markers$pct.2
  markers$gene <- rownames(markers)
  markers$celltype <- cur_celltype
  
  combined_df3 <- rbind(combined_df3, markers)
}
celltype.sALSnoFTLD.markers <- combined_df3

save(celltype.sALSnoFTLD.markers,file="celltype_sALSnoFTLD_markers_LR_RNA.rda")
write.csv(celltype.sALSnoFTLD.markers, file='celltype_sALSnoFTLD_markers_LR_RNA.csv', quote=F)

combined_df4 <- data.frame()
for(i in 1:length(celltypes)){
  
  cur_celltype <- celltypes[[i]]
  print(cur_celltype)
  
  cur_seurat_obj <- subset(intRNA, celltype == cur_celltype)
  Idents(cur_seurat_obj) <- cur_seurat_obj$diagnosis
  print(unique(Idents(cur_seurat_obj)))
  
  markers <- FindMarkers(object = cur_seurat_obj, 
                         ident.1=c('C9ALSFTLD','C9ALSnoFTLD'),
                         ident.2=c('control'),
                         test.use = "LR",
                         latent.vars = c("sex","nCount_RNA"),
                         logfc.threshold = 0.1,
                         min.pct = 0.3,
                         min.diff.pct = -Inf,
                         verbose = TRUE,
                         random.seed = 1,
                         pseudocount.use = 1,
                         assay='RNA')
  markers$diff <- markers$pct.1 - markers$pct.2
  markers$gene <- rownames(markers)
  markers$celltype <- cur_celltype
  
  combined_df4 <- rbind(combined_df4, markers)
}
celltype.C9together.markers <- combined_df3

save(celltype.C9together.markers,file="celltype_C9together_markers_LR_RNA.rda")
write.csv(celltype.C9together.markers, file='celltype_C9together_markers_LR_RNA.csv', quote=F)

####################################
####################################
####################################

library(Seurat)
library(UpSetR)
library(ggplot2)

load("/mnt/WORKHORSE/Jan22_RNA_ATAC/V3only/noRP_DE_logFCpt1/celltype_C9together_markers_MAST_RNA.rda")
load("/mnt/WORKHORSE/Jan22_RNA_ATAC/V3only/noRP_DE_logFCpt1/celltype_C9ALSFTLD_markers_MAST_RNA.rda")
load("/mnt/WORKHORSE/Jan22_RNA_ATAC/V3only/noRP_DE_logFCpt1/celltype_C9ALSnoFTLD_markers_MAST_RNA.rda")
load("/mnt/WORKHORSE/Jan22_RNA_ATAC/V3only/noRP_DE_logFCpt1/celltype_sALSnoFTLD_markers_MAST_RNA.rda")

celltype.C9together.cutoff<-celltype.C9together.markers[celltype.C9together.markers$p_val_adj<0.01, ]
celltype.C9ALSFTLD.cutoff<-celltype.C9ALSFTLD.markers[celltype.C9ALSFTLD.markers$p_val_adj<0.01, ]
celltype.C9ALSnoFTLD.cutoff<-celltype.C9ALSnoFTLD.markers[celltype.C9ALSnoFTLD.markers$p_val_adj<0.01, ]
celltype.sALSnoFTLD.cutoff<-celltype.sALSnoFTLD.markers[celltype.sALSnoFTLD.markers$p_val_adj<0.01, ]

oligo_upset<-list(C9ALSFTLD=celltype.C9ALSFTLD.cutoff[celltype.C9ALSFTLD.cutoff$celltype == "Oligodendrocytes", ]$gene,
                  C9ALSnoFTLD=celltype.C9ALSnoFTLD.cutoff[celltype.C9ALSnoFTLD.cutoff$celltype == "Oligodendrocytes", ]$gene,
                  sALSnoFTLD=celltype.sALSnoFTLD.cutoff[celltype.sALSnoFTLD.cutoff$celltype == "Oligodendrocytes", ]$gene)
#write.csv(astro_upset, file='astro_upset.csv', quote=F)

pdf(file="oligo_upset.pdf", height=4,width=5)
upset(fromList(oligo_upset),group.by="degree",order.by="degree",sets=c("sALSnoFTLD","C9ALSnoFTLD","C9ALSFTLD"),
      sets.bar.color=c("sALSnoFTLD"="#0072B2","C9ALSnoFTLD"="#F0E442","C9ALSFTLD"="#009E73"),
      number.angles=0,point.size=2, line.size=0.5, text.scale=c(1.3,1.3,1.3,1.3,1.3),keep.order=T,
      main.bar.color="gray15", 
      matrix.color=c("gray10"), matrix.dot.alpha=1.5,
      shade.color=c("#0072B2","#009E73"), shade.alpha = 1,
      mainbar.y.label = "Gene Intersection Size",
      sets.x.label = "DE Genes (FDR<0.01)")
dev.off()

opc_upset<-list(C9ALSFTLD=celltype.C9ALSFTLD.cutoff[celltype.C9ALSFTLD.cutoff$celltype == "OPC", ]$gene,
                  C9ALSnoFTLD=celltype.C9ALSnoFTLD.cutoff[celltype.C9ALSnoFTLD.cutoff$celltype == "OPC", ]$gene,
                  sALSnoFTLD=celltype.sALSnoFTLD.cutoff[celltype.sALSnoFTLD.cutoff$celltype == "OPC", ]$gene)
#write.csv(astro_upset, file='astro_upset.csv', quote=F)

pdf(file="opc_upset.pdf", height=4,width=5)
upset(fromList(opc_upset),group.by="degree",order.by="degree",sets=c("sALSnoFTLD","C9ALSnoFTLD","C9ALSFTLD"),
      sets.bar.color=c("sALSnoFTLD"="#0072B2","C9ALSnoFTLD"="#F0E442","C9ALSFTLD"="#009E73"),
      number.angles=0,point.size=2, line.size=0.5, text.scale=c(1.3,1.3,1.3,1.3,1.3),keep.order=T,
      main.bar.color="gray15", 
      matrix.color=c("gray10"), matrix.dot.alpha=1.5,
      shade.color=c("#0072B2","#009E73"), shade.alpha = 1,
      mainbar.y.label = "Gene Intersection Size",
      sets.x.label = "DE Genes (FDR<0.01)")
dev.off()

astro_upset<-list(C9ALSFTLD=celltype.C9ALSFTLD.cutoff[celltype.C9ALSFTLD.cutoff$celltype == "Astrocytes", ]$gene,
                  C9ALSnoFTLD=celltype.C9ALSnoFTLD.cutoff[celltype.C9ALSnoFTLD.cutoff$celltype == "Astrocytes", ]$gene,
                  sALSnoFTLD=celltype.sALSnoFTLD.cutoff[celltype.sALSnoFTLD.cutoff$celltype == "Astrocytes", ]$gene)
#write.csv(astro_upset, file='astro_upset.csv', quote=F)


pdf(file="astro_upset.pdf", height=4,width=5)
upset(fromList(astro_upset),group.by="degree",order.by="degree",sets=c("sALSnoFTLD","C9ALSnoFTLD","C9ALSFTLD"),
      sets.bar.color=c("sALSnoFTLD"="#0072B2","C9ALSnoFTLD"="#F0E442","C9ALSFTLD"="#009E73"),
      number.angles=0,point.size=2, line.size=0.5, text.scale=c(1.3,1.3,1.3,1.3,1.3),keep.order=T,
      main.bar.color="gray15", 
      matrix.color=c("gray10"), matrix.dot.alpha=1.5,
      shade.color=c("#0072B2","#009E73"), shade.alpha = 1,
      mainbar.y.label = "Gene Intersection Size",
      sets.x.label = "DE Genes (FDR<0.01)")
dev.off()

micro_upset<-list(C9ALSFTLD=celltype.C9ALSFTLD.cutoff[celltype.C9ALSFTLD.cutoff$celltype == "Microglia", ]$gene,
               C9ALSnoFTLD=celltype.C9ALSnoFTLD.cutoff[celltype.C9ALSnoFTLD.cutoff$celltype == "Microglia", ]$gene,
               sALSnoFTLD=celltype.sALSnoFTLD.cutoff[celltype.sALSnoFTLD.cutoff$celltype == "Microglia", ]$gene)
#write.csv(astro_upset, file='astro_upset.csv', quote=F)

pdf(file="micro_upset.pdf", height=4,width=5)
upset(fromList(micro_upset),group.by="degree",order.by="degree",sets=c("sALSnoFTLD","C9ALSnoFTLD","C9ALSFTLD"),
      sets.bar.color=c("sALSnoFTLD"="#0072B2","C9ALSnoFTLD"="#F0E442","C9ALSFTLD"="#009E73"),
      number.angles=0,point.size=2, line.size=0.5, text.scale=c(1.3,1.3,1.3,1.3,1.3),keep.order=T,
      main.bar.color="gray15", 
      matrix.color=c("gray10"), matrix.dot.alpha=1.5,
      shade.color=c("#0072B2","#009E73"), shade.alpha = 1,
      mainbar.y.label = "Gene Intersection Size",
      sets.x.label = "DE Genes (FDR<0.01)")
dev.off()

ex_upset<-list(C9ALSFTLD=celltype.C9ALSFTLD.cutoff[celltype.C9ALSFTLD.cutoff$celltype == "Excitatory", ]$gene,
                  C9ALSnoFTLD=celltype.C9ALSnoFTLD.cutoff[celltype.C9ALSnoFTLD.cutoff$celltype == "Excitatory", ]$gene,
                  sALSnoFTLD=celltype.sALSnoFTLD.cutoff[celltype.sALSnoFTLD.cutoff$celltype == "Excitatory", ]$gene)
#write.csv(astro_upset, file='astro_upset.csv', quote=F)

pdf(file="ex_upset.pdf", height=4,width=5)
upset(fromList(ex_upset),group.by="degree",order.by="degree",sets=c("sALSnoFTLD","C9ALSnoFTLD","C9ALSFTLD"),
      sets.bar.color=c("sALSnoFTLD"="#0072B2","C9ALSnoFTLD"="#F0E442","C9ALSFTLD"="#009E73"),
      number.angles=0,point.size=2, line.size=0.5, text.scale=c(1.3,1.3,1.3,1.3,1.3),keep.order=T,
      main.bar.color="gray15", 
      matrix.color=c("gray10"), matrix.dot.alpha=1.5,
      shade.color=c("#0072B2","#009E73"), shade.alpha = 1,
      mainbar.y.label = "Gene Intersection Size",
      sets.x.label = "DE Genes (FDR<0.01)")
dev.off()

in_upset<-list(C9ALSFTLD=celltype.C9ALSFTLD.cutoff[celltype.C9ALSFTLD.cutoff$celltype == "Inhibitory", ]$gene,
               C9ALSnoFTLD=celltype.C9ALSnoFTLD.cutoff[celltype.C9ALSnoFTLD.cutoff$celltype == "Inhibitory", ]$gene,
               sALSnoFTLD=celltype.sALSnoFTLD.cutoff[celltype.sALSnoFTLD.cutoff$celltype == "Inhibitory", ]$gene)
#write.csv(astro_upset, file='astro_upset.csv', quote=F)

pdf(file="in_upset.pdf", height=4,width=5)
upset(fromList(in_upset),group.by="degree",order.by="degree",sets=c("sALSnoFTLD","C9ALSnoFTLD","C9ALSFTLD"),
      sets.bar.color=c("sALSnoFTLD"="#0072B2","C9ALSnoFTLD"="#F0E442","C9ALSFTLD"="#009E73"),
      number.angles=0,point.size=2, line.size=0.5, text.scale=c(1.3,1.3,1.3,1.3,1.3),keep.order=T,
      main.bar.color="gray15", 
      matrix.color=c("gray10"), matrix.dot.alpha=1.5,
      shade.color=c("#0072B2","#009E73"), shade.alpha = 1,
      mainbar.y.label = "Gene Intersection Size",
      sets.x.label = "DE Genes (FDR<0.01)")
dev.off()

#####################
## by cell subtype
#####################

load("/mnt/WORKHORSE/Jan22_RNA_ATAC/V3only/noRP_DE_logFCpt1/cellsubtype_C9together_markers_LR_RNA.rda")
load("/mnt/WORKHORSE/Jan22_RNA_ATAC/V3only/noRP_DE_logFCpt1/cellsubtype_C9ALSFTLD_markers_LR_RNA.rda")
load("/mnt/WORKHORSE/Jan22_RNA_ATAC/V3only/noRP_DE_logFCpt1/cellsubtype_C9ALSnoFTLD_markers_LR_RNA.rda")
load("/mnt/WORKHORSE/Jan22_RNA_ATAC/V3only/noRP_DE_logFCpt1/cellsubtype_sALSnoFTLD_markers_LR_RNA.rda")

cellsubtype.C9together.cutoff<-cellsubtype.C9together.markers[cellsubtype.C9together.markers$p_val_adj<0.01, ]
cellsubtype.C9ALSFTLD.cutoff<-cellsubtype.C9ALSFTLD.markers[cellsubtype.C9ALSFTLD.markers$p_val_adj<0.01, ]
cellsubtype.C9ALSnoFTLD.cutoff<-cellsubtype.C9ALSnoFTLD.markers[cellsubtype.C9ALSnoFTLD.markers$p_val_adj<0.01, ]
cellsubtype.sALSnoFTLD.cutoff<-cellsubtype.sALSnoFTLD.markers[cellsubtype.sALSnoFTLD.markers$p_val_adj<0.01, ]


ASTPP_upset<-list(C9ALSFTLD=cellsubtype.C9ALSFTLD.cutoff[cellsubtype.C9ALSFTLD.cutoff$cellsubtype == "AST-PP", ]$gene,
                  C9ALSnoFTLD=cellsubtype.C9ALSnoFTLD.cutoff[cellsubtype.C9ALSnoFTLD.cutoff$cellsubtype == "AST-PP", ]$gene,
                  sALSnoFTLD=cellsubtype.sALSnoFTLD.cutoff[cellsubtype.sALSnoFTLD.cutoff$cellsubtype == "AST-PP", ]$gene)
#write.csv(astro_upset, file='astro_upset.csv', quote=F)

pdf(file="ASTPP_upset.pdf", height=4,width=5)
upset(fromList(ASTPP_upset),group.by="degree",order.by="degree",sets=c("sALSnoFTLD","C9ALSnoFTLD","C9ALSFTLD"),
      sets.bar.color=c("sALSnoFTLD"="#0072B2","C9ALSnoFTLD"="#F0E442","C9ALSFTLD"="#009E73"),
      number.angles=0,point.size=2, line.size=0.5, text.scale=c(1.3,1.3,1.3,1.3,1.3),keep.order=T,
      main.bar.color="gray15", 
      matrix.color=c("gray10"), matrix.dot.alpha=1.5,
      shade.color=c("#0072B2","#009E73"), shade.alpha = 1,
      mainbar.y.label = "Gene Intersection Size",
      sets.x.label = "DE Genes (FDR<0.01)")
dev.off()

ASTFB_upset<-list(C9ALSFTLD=cellsubtype.C9ALSFTLD.cutoff[cellsubtype.C9ALSFTLD.cutoff$cellsubtype == "AST-FB", ]$gene,
                  C9ALSnoFTLD=cellsubtype.C9ALSnoFTLD.cutoff[cellsubtype.C9ALSnoFTLD.cutoff$cellsubtype == "AST-FB", ]$gene,
                  sALSnoFTLD=cellsubtype.sALSnoFTLD.cutoff[cellsubtype.sALSnoFTLD.cutoff$cellsubtype == "AST-FB", ]$gene)
#write.csv(astro_upset, file='astro_upset.csv', quote=F)

pdf(file="ASTFB_upset.pdf", height=4,width=5)
upset(fromList(ASTFB_upset),group.by="degree",order.by="degree",sets=c("sALSnoFTLD","C9ALSnoFTLD","C9ALSFTLD"),
      sets.bar.color=c("sALSnoFTLD"="#0072B2","C9ALSnoFTLD"="#F0E442","C9ALSFTLD"="#009E73"),
      number.angles=0,point.size=2, line.size=0.5, text.scale=c(1.3,1.3,1.3,1.3,1.3),keep.order=T,
      main.bar.color="gray15", 
      matrix.color=c("gray10"), matrix.dot.alpha=1.5,
      shade.color=c("#0072B2","#009E73"), shade.alpha = 1,
      mainbar.y.label = "Gene Intersection Size",
      sets.x.label = "DE Genes (FDR<0.01)")
dev.off()

L23_upset<-list(C9ALSFTLD=cellsubtype.C9ALSFTLD.cutoff[cellsubtype.C9ALSFTLD.cutoff$cellsubtype == "L2/3", ]$gene,
                  C9ALSnoFTLD=cellsubtype.C9ALSnoFTLD.cutoff[cellsubtype.C9ALSnoFTLD.cutoff$cellsubtype == "L2/3", ]$gene,
                  sALSnoFTLD=cellsubtype.sALSnoFTLD.cutoff[cellsubtype.sALSnoFTLD.cutoff$cellsubtype == "L2/3", ]$gene)
#write.csv(astro_upset, file='astro_upset.csv', quote=F)

pdf(file="L23_upset.pdf", height=4,width=5)
upset(fromList(L23_upset),group.by="degree",order.by="degree",sets=c("sALSnoFTLD","C9ALSnoFTLD","C9ALSFTLD"),
      sets.bar.color=c("sALSnoFTLD"="#0072B2","C9ALSnoFTLD"="#F0E442","C9ALSFTLD"="#009E73"),
      number.angles=0,point.size=2, line.size=0.5, text.scale=c(1.3,1.3,1.3,1.3,1.3),keep.order=T,
      main.bar.color="gray15", 
      matrix.color=c("gray10"), matrix.dot.alpha=1.5,
      shade.color=c("#0072B2","#009E73"), shade.alpha = 1,
      mainbar.y.label = "Gene Intersection Size",
      sets.x.label = "DE Genes (FDR<0.01)")
dev.off()

L4_upset<-list(C9ALSFTLD=cellsubtype.C9ALSFTLD.cutoff[cellsubtype.C9ALSFTLD.cutoff$cellsubtype == "L4", ]$gene,
                  C9ALSnoFTLD=cellsubtype.C9ALSnoFTLD.cutoff[cellsubtype.C9ALSnoFTLD.cutoff$cellsubtype == "L4", ]$gene,
                  sALSnoFTLD=cellsubtype.sALSnoFTLD.cutoff[cellsubtype.sALSnoFTLD.cutoff$cellsubtype == "L4", ]$gene)
#write.csv(astro_upset, file='astro_upset.csv', quote=F)

pdf(file="L4_upset.pdf", height=4,width=5)
upset(fromList(L4_upset),group.by="degree",order.by="degree",sets=c("sALSnoFTLD","C9ALSnoFTLD","C9ALSFTLD"),
      sets.bar.color=c("sALSnoFTLD"="#0072B2","C9ALSnoFTLD"="#F0E442","C9ALSFTLD"="#009E73"),
      number.angles=0,point.size=2, line.size=0.5, text.scale=c(1.3,1.3,1.3,1.3,1.3),keep.order=T,
      main.bar.color="gray15", 
      matrix.color=c("gray10"), matrix.dot.alpha=1.5,
      shade.color=c("#0072B2","#009E73"), shade.alpha = 1,
      mainbar.y.label = "Gene Intersection Size",
      sets.x.label = "DE Genes (FDR<0.01)")
dev.off()

L56_upset<-list(C9ALSFTLD=cellsubtype.C9ALSFTLD.cutoff[cellsubtype.C9ALSFTLD.cutoff$cellsubtype == "L5/6", ]$gene,
                  C9ALSnoFTLD=cellsubtype.C9ALSnoFTLD.cutoff[cellsubtype.C9ALSnoFTLD.cutoff$cellsubtype == "L5/6", ]$gene,
                  sALSnoFTLD=cellsubtype.sALSnoFTLD.cutoff[cellsubtype.sALSnoFTLD.cutoff$cellsubtype == "L5/6", ]$gene)
#write.csv(astro_upset, file='astro_upset.csv', quote=F)

pdf(file="L56_upset.pdf", height=4,width=5)
upset(fromList(L56_upset),group.by="degree",order.by="degree",sets=c("sALSnoFTLD","C9ALSnoFTLD","C9ALSFTLD"),
      sets.bar.color=c("sALSnoFTLD"="#0072B2","C9ALSnoFTLD"="#F0E442","C9ALSFTLD"="#009E73"),
      number.angles=0,point.size=2, line.size=0.5, text.scale=c(1.3,1.3,1.3,1.3,1.3),keep.order=T,
      main.bar.color="gray15", 
      matrix.color=c("gray10"), matrix.dot.alpha=1.5,
      shade.color=c("#0072B2","#009E73"), shade.alpha = 1,
      mainbar.y.label = "Gene Intersection Size",
      sets.x.label = "DE Genes (FDR<0.01)")
dev.off()

L56CC_upset<-list(C9ALSFTLD=cellsubtype.C9ALSFTLD.cutoff[cellsubtype.C9ALSFTLD.cutoff$cellsubtype == "L5/6-CC", ]$gene,
                  C9ALSnoFTLD=cellsubtype.C9ALSnoFTLD.cutoff[cellsubtype.C9ALSnoFTLD.cutoff$cellsubtype == "L5/6-CC", ]$gene,
                  sALSnoFTLD=cellsubtype.sALSnoFTLD.cutoff[cellsubtype.sALSnoFTLD.cutoff$cellsubtype == "L5/6-CC", ]$gene)
#write.csv(astro_upset, file='astro_upset.csv', quote=F)

pdf(file="L56cc_upset.pdf", height=4,width=5)
upset(fromList(L56CC_upset),group.by="degree",order.by="degree",sets=c("sALSnoFTLD","C9ALSnoFTLD","C9ALSFTLD"),
      sets.bar.color=c("sALSnoFTLD"="#0072B2","C9ALSnoFTLD"="#F0E442","C9ALSFTLD"="#009E73"),
      number.angles=0,point.size=2, line.size=0.5, text.scale=c(1.3,1.3,1.3,1.3,1.3),keep.order=T,
      main.bar.color="gray15", 
      matrix.color=c("gray10"), matrix.dot.alpha=1.5,
      shade.color=c("#0072B2","#009E73"), shade.alpha = 1,
      mainbar.y.label = "Gene Intersection Size",
      sets.x.label = "DE Genes (FDR<0.01)")
dev.off()

INPV_upset<-list(C9ALSFTLD=cellsubtype.C9ALSFTLD.cutoff[cellsubtype.C9ALSFTLD.cutoff$cellsubtype == "IN-PV", ]$gene,
                  C9ALSnoFTLD=cellsubtype.C9ALSnoFTLD.cutoff[cellsubtype.C9ALSnoFTLD.cutoff$cellsubtype == "IN-PV", ]$gene,
                  sALSnoFTLD=cellsubtype.sALSnoFTLD.cutoff[cellsubtype.sALSnoFTLD.cutoff$cellsubtype == "IN-PV", ]$gene)
#write.csv(astro_upset, file='astro_upset.csv', quote=F)

pdf(file="INPV_upset.pdf", height=4,width=5)
upset(fromList(INPV_upset),group.by="degree",order.by="degree",sets=c("sALSnoFTLD","C9ALSnoFTLD","C9ALSFTLD"),
      sets.bar.color=c("sALSnoFTLD"="#0072B2","C9ALSnoFTLD"="#F0E442","C9ALSFTLD"="#009E73"),
      number.angles=0,point.size=2, line.size=0.5, text.scale=c(1.3,1.3,1.3,1.3,1.3),keep.order=T,
      main.bar.color="gray15", 
      matrix.color=c("gray10"), matrix.dot.alpha=1.5,
      shade.color=c("#0072B2","#009E73"), shade.alpha = 1,
      mainbar.y.label = "Gene Intersection Size",
      sets.x.label = "DE Genes (FDR<0.01)")
dev.off()

INSST_upset<-list(C9ALSFTLD=cellsubtype.C9ALSFTLD.cutoff[cellsubtype.C9ALSFTLD.cutoff$cellsubtype == "IN-SST", ]$gene,
                  C9ALSnoFTLD=cellsubtype.C9ALSnoFTLD.cutoff[cellsubtype.C9ALSnoFTLD.cutoff$cellsubtype == "IN-SST", ]$gene,
                  sALSnoFTLD=cellsubtype.sALSnoFTLD.cutoff[cellsubtype.sALSnoFTLD.cutoff$cellsubtype == "IN-SST", ]$gene)
#write.csv(astro_upset, file='astro_upset.csv', quote=F)

pdf(file="INSST_upset.pdf", height=4,width=5)
upset(fromList(INSST_upset),group.by="degree",order.by="degree",sets=c("sALSnoFTLD","C9ALSnoFTLD","C9ALSFTLD"),
      sets.bar.color=c("sALSnoFTLD"="#0072B2","C9ALSnoFTLD"="#F0E442","C9ALSFTLD"="#009E73"),
      number.angles=0,point.size=2, line.size=0.5, text.scale=c(1.3,1.3,1.3,1.3,1.3),keep.order=T,
      main.bar.color="gray15", 
      matrix.color=c("gray10"), matrix.dot.alpha=1.5,
      shade.color=c("#0072B2","#009E73"), shade.alpha = 1,
      mainbar.y.label = "Gene Intersection Size",
      sets.x.label = "DE Genes (FDR<0.01)")
dev.off()

INSV2C_upset<-list(C9ALSFTLD=cellsubtype.C9ALSFTLD.cutoff[cellsubtype.C9ALSFTLD.cutoff$cellsubtype == "IN-SV2C", ]$gene,
                  C9ALSnoFTLD=cellsubtype.C9ALSnoFTLD.cutoff[cellsubtype.C9ALSnoFTLD.cutoff$cellsubtype == "IN-SV2C", ]$gene,
                  sALSnoFTLD=cellsubtype.sALSnoFTLD.cutoff[cellsubtype.sALSnoFTLD.cutoff$cellsubtype == "IN-SV2C", ]$gene)
#write.csv(astro_upset, file='astro_upset.csv', quote=F)

pdf(file="INSV2C_upset.pdf", height=4,width=5)
upset(fromList(INSV2C_upset),group.by="degree",order.by="degree",sets=c("sALSnoFTLD","C9ALSnoFTLD","C9ALSFTLD"),
      sets.bar.color=c("sALSnoFTLD"="#0072B2","C9ALSnoFTLD"="#F0E442","C9ALSFTLD"="#009E73"),
      number.angles=0,point.size=2, line.size=0.5, text.scale=c(1.3,1.3,1.3,1.3,1.3),keep.order=T,
      main.bar.color="gray15", 
      matrix.color=c("gray10"), matrix.dot.alpha=1.5,
      shade.color=c("#0072B2","#009E73"), shade.alpha = 1,
      mainbar.y.label = "Gene Intersection Size",
      sets.x.label = "DE Genes (FDR<0.01)")

INVIP_upset<-list(C9ALSFTLD=cellsubtype.C9ALSFTLD.cutoff[cellsubtype.C9ALSFTLD.cutoff$cellsubtype == "IN-VIP", ]$gene,
                  C9ALSnoFTLD=cellsubtype.C9ALSnoFTLD.cutoff[cellsubtype.C9ALSnoFTLD.cutoff$cellsubtype == "IN-VIP", ]$gene,
                  sALSnoFTLD=cellsubtype.sALSnoFTLD.cutoff[cellsubtype.sALSnoFTLD.cutoff$cellsubtype == "IN-VIP", ]$gene)
#write.csv(astro_upset, file='astro_upset.csv', quote=F)

pdf(file="INVIP_upset.pdf", height=4,width=5)
upset(fromList(INVIP_upset),group.by="degree",order.by="degree",sets=c("sALSnoFTLD","C9ALSnoFTLD","C9ALSFTLD"),
      sets.bar.color=c("sALSnoFTLD"="#0072B2","C9ALSnoFTLD"="#F0E442","C9ALSFTLD"="#009E73"),
      number.angles=0,point.size=2, line.size=0.5, text.scale=c(1.3,1.3,1.3,1.3,1.3),keep.order=T,
      main.bar.color="gray15", 
      matrix.color=c("gray10"), matrix.dot.alpha=1.5,
      shade.color=c("#0072B2","#009E73"), shade.alpha = 1,
      mainbar.y.label = "Gene Intersection Size",
      sets.x.label = "DE Genes (FDR<0.01)")
dev.off()

library(fgsea)
library(gmt)


goALL_gmt <- gmtPathways("genesets/Human_GOBP_AllPathways_no_GO_iea_January_01_2022_symbol.gmt")
goBP_gmt <- gmtPathways("genesets/Human_GO_bp_no_GO_iea_symbol.gmt")
goBP_gmt

ranks_in<-cellsubtype.C9ALSFTLD.cutoff[cellsubtype.C9ALSFTLD.cutoff$cellsubtype == "IN-VIP", ]$avg_log2FC
names(ranks_in)<-cellsubtype.C9ALSFTLD.cutoff[cellsubtype.C9ALSFTLD.cutoff$cellsubtype == "IN-VIP", ]$gene

fgseaRes_in<-fgsea(pathways=goBP_gmt,stats=ranks_in,minSize=10,maxSize=1000)
topPathwaysUp <- fgseaRes_in[ES > 0][head(order(pval), n=20), pathway]
topPathwaysDown <- fgseaRes_in[ES < 0][head(order(pval), n=20), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
topPathways

##########################
### subsetting based on celltype:
##########################


Idents(intRNA)<-"celltype"
intRNA.astro<-subset(intRNA,idents="Astrocytes")
intRNA.microglia<-subset(intRNA,idents="Microglia")
intRNA.opc<-subset(intRNA,idents="OPC")
intRNA.oligo<-subset(intRNA,idents="Oligodendrocytes")
intRNA.excitatory<-subset(intRNA,idents="Excitatory")
intRNA.inhibitory<-subset(intRNA,idents="Inhibitory")
Idents(intRNA)<-"cellsubtype"
intRNA.astro1<-subset(intRNA,idents="AST-PP")
intRNA.astro2<-subset(intRNA,idents="AST-FB")
intRNA.ex1<-subset(intRNA,idents="L23")
intRNA.ex2<-subset(intRNA,idents="L4")
intRNA.ex3<-subset(intRNA,idents="L56")
intRNA.ex4<-subset(intRNA,idents="L56-CC")
intRNA.in1<-subset(intRNA,idents="IN-PV")
intRNA.in2<-subset(intRNA,idents="IN-SST")
intRNA.in3<-subset(intRNA,idents="IN-SV2C")
intRNA.in4<-subset(intRNA,idents="IN-VIP")

#save files
saveRDS(intRNA.excitatory, file="objects/V3only_noRP_excitatory.RDS")
saveRDS(intRNA.inhibitory, file="objects/V3only_noRP_inhibitory.RDS")
saveRDS(intRNA.microglia, file="objects/V3only_noRP_micro.RDS")
saveRDS(intRNA.opc, file="objects/V3only_noRP_OPC.RDS")
saveRDS(intRNA.oligo, file="objects/V3only_noRP_oligo.RDS")
saveRDS(intRNA.astro, file="objects/V3only_noRP_astro.RDS")
saveRDS(intRNA.astro1, file="objects/V3only_noRP_AST-PP.RDS")
saveRDS(intRNA.astro2, file="objects/V3only_noRP_AST-FB.RDS")
saveRDS(intRNA.ex1,file="objects/V3only_noRP_L2-3.RDS")
saveRDS(intRNA.ex2,file="objects/V3only_noRP_L4.RDS")
saveRDS(intRNA.ex3,file="objects/V3only_noRP_L5-6.RDS")
saveRDS(intRNA.ex4,file="objects/V3only_noRP_L5-6CC.RDS")
saveRDS(intRNA.in1,file="objects/V3only_noRP_IN-PV.RDS")
saveRDS(intRNA.in2,file="objects/V3only_noRP_IN-SST.RDS")
saveRDS(intRNA.in3,file="objects/V3only_noRP_IN-SV2C.RDS")
saveRDS(intRNA.in4,file="objects/V3only_noRP_IN-VIP.RDS")