library(Seurat)
set.seed(1234)

intRNA<-readRDS(file="objects/intRNA_labelxfer.RDS")

intRNA$sample<-factor(intRNA$sample, levels=c("CTRL1","CTRL2","CTRL3","CTRL4","CTRL5","CTRL6","C9noALSnoFTLD",
                                              "C9ALSFTLD1","C9ALSFTLD2","C9ALSFTLD3","C9ALSFTLD4","C9ALSFTLD5","C9ALSFTLD6","C9ALSnoFTLD1","C9ALSnoFTLD2","C9ALSnoFTLD3",
                                              "sALSnoFTLD1","sALSnoFTLD2","sALSnoFTLD3","sALSnoFTLD4","sALSnoFTLD5","sALSnoFTLD6","sALSnoFTLD7","sALSnoFTLD8"))

intRNA$cellsubtype <- factor(intRNA$cellsubtype, levels = c("Oligodendrocytes","OPC","AST-PP","AST-FB","Microglia","Endothelial",
                                                                      "L2/3","L4","L5/6","L5/6-CC","Neu-NRGN-I","Neu-NRGN-II","Neu-mat","IN-PV","IN-SST","IN-SV2C","IN-VIP"))

#controls PFC Velmeshev et al
intPFC$cluster <- factor(intPFC$cluster, levels = c("Oligodendrocytes","OPC","AST-PP","AST-FB","Microglia","Endothelial",
                                                                      "L2/3","L4","L5/6","L5/6-CC","Neu-NRGN-I","Neu-NRGN-II","Neu-mat","IN-PV","IN-SST","IN-SV2C","IN-VIP"))

intRNA$diagnosis <- factor(intRNA$diagnosis, levels = c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))

intRNA$diagnosis_cellsubtype <- paste0(intRNA$diagnosis, "_", intRNA$cellsubtype)

Idents(intRNA)<-"cellsubtype"
celltypes<-c("Oligodendrocytes","OPC","Astrocytes","Astrocytes","Microglia","Endothelial",
"Excitatory","Excitatory","Excitatory","Excitatory","Excitatory",
"Inhibitory","Inhibitory","Inhibitory","Inhibitory")
names(celltypes)<-levels(intRNA)
intRNA<-RenameIdents(intRNA,celltypes)
intRNA$celltype<-Idents(intRNA)
intRNA$diagnosis_celltype <- paste0(intRNA$diagnosis, "_", intRNA$celltype)

intRNA <- RenameIdents(intRNA, `L2/3` = "L23", `L5/6` = "L56", `L5/6-CC` = "L56-CC")
Idents(intRNA) <- factor(Idents(intRNA), levels = c("Oligodendrocytes","OPC","AST-PP","AST-FB","Microglia",
                                                    "L23","L4","L56","L56-CC","IN-PV","IN-SST","IN-SV2C","IN-VIP"))
intRNA$cellsubtype<-Idents(intRNA)

intRNA$diagnosis_celltype<-factor(intRNA$diagnosis_celltype,c(
  "control_Oligodendrocytes","C9noALSnoFTLD_Oligodendrocytes","C9ALSFTLD_Oligodendrocytes","C9ALSnoFTLD_Oligodendrocytes","sALSnoFTLD_Oligodendrocytes",
  "control_OPC","C9noALSnoFTLD_OPC","C9ALSFTLD_OPC","C9ALSnoFTLD_OPC","sALSnoFTLD_OPC",
  "control_Astrocytes","C9noALSnoFTLD_Astrocytes","C9ALSFTLD_Astrocytes","C9ALSnoFTLD_Astrocytes","sALSnoFTLD_Astrocytes",
  "control_Microglia","C9noALSnoFTLD_Microglia","C9ALSFTLD_Microglia","C9ALSnoFTLD_Microglia","sALSnoFTLD_Microglia",
  "control_Excitatory","C9noALSnoFTLD_Excitatory","C9ALSFTLD_Excitatory","C9ALSnoFTLD_Excitatory","sALSnoFTLD_Excitatory",
  "control_Inhibitory","C9noALSnoFTLD_Inhibitory","C9ALSFTLD_Inhibitory","C9ALSnoFTLD_Inhibitory","sALSnoFTLD_Inhibitory"
  ))

intRNA$diagnosis_cellsubtype<-factor(intRNA$diagnosis_cellsubtype,c(
  "control_Oligodendrocytes","C9noALSnoFTLD_Oligodendrocytes","C9ALSFTLD_Oligodendrocytes","C9ALSnoFTLD_Oligodendrocytes","sALSnoFTLD_Oligodendrocytes",
  "control_OPC","C9noALSnoFTLD_OPC","C9ALSFTLD_OPC","C9ALSnoFTLD_OPC","sALSnoFTLD_OPC",
  "control_AST-FB","C9noALSnoFTLD_AST-FB","C9ALSFTLD_AST-FB","C9ALSnoFTLD_AST-FB","sALSnoFTLD_AST-FB",
  "control_AST-PP","C9noALSnoFTLD_AST-PP","C9ALSFTLD_AST-PP","C9ALSnoFTLD_AST-PP","sALSnoFTLD_AST-PP",
  "control_Microglia","C9noALSnoFTLD_Microglia","C9ALSFTLD_Microglia","C9ALSnoFTLD_Microglia","sALSnoFTLD_Microglia",
  "control_L23","C9noALSnoFTLD_L23","C9ALSFTLD_L23","C9ALSnoFTLD_L23","sALSnoFTLD_L23",
  "control_L4","C9noALSnoFTLD_L4","C9ALSFTLD_L4","C9ALSnoFTLD_L4","sALSnoFTLD_L4",
  "control_L56","C9noALSnoFTLD_L56","C9ALSFTLD_L56","C9ALSnoFTLD_L56","sALSnoFTLD_L56",
  "control_L56-CC","C9noALSnoFTLD_L56-CC","C9ALSFTLD_L56-CC","C9ALSnoFTLD_L56-CC","sALSnoFTLD_L56-CC",
  "control_IN-PV","C9noALSnoFTLD_IN-PV","C9ALSFTLD_IN-PV","C9ALSnoFTLD_IN-PV","sALSnoFTLD_IN-PV",
  "control_IN-SST","C9noALSnoFTLD_IN-SST","C9ALSFTLD_IN-SST","C9ALSnoFTLD_IN-SST","sALSnoFTLD_IN-SST",
  "control_IN-SV2C","C9noALSnoFTLD_IN-SV2C","C9ALSFTLD_IN-SV2C","C9ALSnoFTLD_IN-SV2C","sALSnoFTLD_IN-SV2C",
  "control_IN-VIP","C9noALSnoFTLD_IN-VIP","C9ALSFTLD_IN-VIP","C9ALSnoFTLD_IN-VIP","sALSnoFTLD_IN-VIP"
  ))

saveRDS(intRNA,file="objects/intRNA_final.RDS")

##remove ribo genes for DE analysis based on over-enrichment when included in dataset
counts<-GetAssayData(intRNA,assay="RNA")
counts
ribo.genes <- c(rownames(intRNA@assays$RNA)[grep("^RP[1:9]", rownames(intRNA@assays$RNA))],
                rownames(intRNA@assays$RNA)[grep("^RP[L,S]", rownames(intRNA@assays$RNA))]
)
ribo.genes
counts<-counts[-(which(rownames(counts) %in% ribo.genes)),]
intRNA2<-subset(intRNA,features=rownames(counts))
intRNA2
rm(intRNA)
gc()
intRNA<-intRNA2
rm(intRNA2)
gc()

Idents(intRNA)<-"diagnosis"
Idents(intRNA)<-"pseudobulk"
pseudobulk<-c("pseudobulk","pseudobulk","pseudobulk","pseudobulk","pseudobulk")
names(pseudobulk)<-levels(intRNA)
intRNA<-RenameIdents(intRNA,pseudobulk)
intRNA$pseudobulk<-Idents(intRNA)

pseudobulk<-c("snRNA")
names(pseudobulk)<-levels(intRNA)
intRNA<-RenameIdents(intRNA,pseudobulk)
intRNA$pseudobulk<-Idents(intRNA)

saveRDS(intRNA, file="objects/intRNA_final.RDS")


# subset by V3 chemistry
Idents(intRNA)<-"sample"
intRNA2<-subset(intRNA,idents=c('CTRL1','CTRL2','CTRL3','CTRL4','CTRL5','CTRL6','C9noALSnoFTLD','C9ALSnoFTLD1','C9ALSnoFTLD3','C9ALSFTLD1','C9ALSFTLD2','C9ALSFTLD3','sALSnoFTLD3','sALSnoFTLD4','sALSnoFTLD5','sALSnoFTLD8'))
intRNA2$sample<-factor(intRNA2$sample,levels=c('CTRL1','CTRL2','CTRL3','CTRL4','CTRL5','CTRL6','C9noALSnoFTLD','C9ALSnoFTLD1','C9ALSnoFTLD3','C9ALSFTLD1','C9ALSFTLD2','C9ALSFTLD3','sALSnoFTLD3','sALSnoFTLD4','sALSnoFTLD5','sALSnoFTLD8'))
table(intRNA2$sample)
saveRDS(intRNA2,file="intRNA_V3only.RDS")

DefaultAssay(intRNA)<-"RNA"
Idents(intRNA)<-"celltype"
table(intRNA$celltype)
intRNA.oligo<-subset(intRNA, idents=c("Oligodendrocytes"))
saveRDS(intRNA.oligo, file="objects/RNA_oligo_noRP.RDS")

intRNA.opc<-subset(intRNA, idents=c("OPC"))
saveRDS(intRNA.opc, file="objects/RNA_OPCs_noRP.RDS")

intRNA.astro<-subset(intRNA, idents=c("Astrocytes"))
saveRDS(intRNA.astro, file="objects/RNA_astrocytes_noRP.RDS")

intRNA.microglia<-subset(intRNA, idents=c("Microglia"))
saveRDS(intRNA.microglia, file="objects/RNA_microglia_noRP.RDS")

intRNA.endo<-subset(intRNA, idents=c("Endothelial"))
saveRDS(intRNA.endo, file="objects/RNA_endo_noRP.RDS")

intRNA.excitatory<-subset(intRNA, idents=c("Excitatory"))
saveRDS(intRNA.excitatory, file="objects/RNA_excitatory_noRP.RDS")

intRNA.inhibitory<-subset(intRNA, idents=c("Inhibitory"))
saveRDS(intRNA.inhibitory, file="objects/RNA_inhibitory_noRP.RDS")

intRNA

Idents(intRNA)<-"cellsubtype"
intRNA.astro1<-subset(intRNA, idents=c("AST-PP"))
intRNA.astro2<-subset(intRNA, idents=c("AST-FB"))
saveRDS(intRNA.astro1, file="objects/RNA_AST-PP_noRP.RDS")
saveRDS(intRNA.astro2, file="objects/RNA_AST-FB_noRP.RDS")

Idents(intRNA)<-"cellsubtype"
intRNA.ex1<-subset(intRNA, idents=c("L2/3"))
intRNA.ex2<-subset(intRNA, idents=c("L4"))
intRNA.ex3<-subset(intRNA, idents=c("L5/6"))
intRNA.ex4<-subset(intRNA, idents=c("L5/6-CC"))
saveRDS(intRNA.ex1, file="objects/RNA_L2-3_noRP.RDS")
saveRDS(intRNA.ex2, file="objects/RNA_L4_noRP.RDS")
saveRDS(intRNA.ex3, file="objects/RNA_L5-6_noRP.RDS")
saveRDS(intRNA.ex4, file="objects/RNA_L5-6CC_noRP.RDS")


Idents(intRNA)<-"cellsubtype"
intRNA.in1<-subset(intRNA, idents=c("IN-PV"))
intRNA.in2<-subset(intRNA, idents=c("IN-SST"))
intRNA.in3<-subset(intRNA, idents=c("IN-SV2C"))
intRNA.in4<-subset(intRNA, idents=c("IN-VIP"))
saveRDS(intRNA.in1, file="objects/RNA_IN-PV_noRP.RDS")
saveRDS(intRNA.in2, file="objects/RNA_IN-SST_noRP.RDS")
saveRDS(intRNA.in3, file="objects/RNA_IN-SV2C_noRP.RDS")
saveRDS(intRNA.in4, file="objects/RNA_IN-VIP_noRP.RDS")

