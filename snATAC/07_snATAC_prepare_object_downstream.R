library(Seurat)
library(Signac)
set.seed(1234)

intATAC<-readRDS(file="objects/ATAC_SignacPeaks.RDS")

fragments<-CreateFragmentObject(path="fragments/fragments.tsv.gz",cells=colnames(intATAC))
Fragments(intATAC)<-NULL
Fragments(intATAC)<-fragments
DefaultAssay(intATAC) <- "ATAC"

intATAC$diagnosis<-factor(intATAC$diagnosis,levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))

Idents(intATAC)<-"cellsubtype"


Idents(intATAC)<-"cellsubtype"
celltypes<-c("Oligodendrocytes","OPC","Astrocytes","Astrocytes",
    "Microglia","Endothelial","Excitatory","Excitatory","Excitatory",
    "Excitatory","Inhibitory","Inhibitory","Inhibitory","Inhibitory")
names(celltypes)<-levels(intATAC)
intATAC<-RenameIdents(intATAC,celltypes)
intATAC$celltype<-Idents(intATAC)
intATAC$diagnosis_celltype <- paste0(intATAC$diagnosis, "_", intATAC$celltype)

intATAC<-subset(intATAC,idents=c("Oligodendrocytes","OPC","AST-PP","AST-FB","Microglia",
                                 "L2/3","L4","L5/6","L5/6","IN-PV","IN-SST","IN-SV2C","IN-VIP"))
Idents(intATAC)<-"diagnosis"

intATAC$diagnosis<-factor(intATAC$diagnosis,levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))

intRNA$diagnosis_celltype <- paste0(intRNA$diagnosis, "_", intRNA$celltype)

intATAC <- RenameIdents(intATAC, `L2/3` = "L23", `L5/6` = "L56", `L5/6-CC` = "L56-CC")
Idents(intATAC) <- factor(Idents(intATAC), levels = c("Oligodendrocytes","OPC","AST-PP","AST-FB","Microglia",
                                                    "L23","L4","L56","L56-CC","IN-PV","IN-SST","IN-SV2C","IN-VIP"))
intATAC$cellsubtype<-Idents(intATAC)

intATAC$diagnosis_celltype <- paste0(intATAC$diagnosis, "_", intATAC$celltype)

intATAC$diagnosis_celltype<-factor(intATAC$diagnosis_celltype,c(
  "control_Oligodendrocytes","C9noALSnoFTLD_Oligodendrocytes","C9ALSFTLD_Oligodendrocytes","C9ALSnoFTLD_Oligodendrocytes","sALSnoFTLD_Oligodendrocytes",
  "control_OPC","C9noALSnoFTLD_OPC","C9ALSFTLD_OPC","C9ALSnoFTLD_OPC","sALSnoFTLD_OPC",
  "control_Astrocytes","C9noALSnoFTLD_Astrocytes","C9ALSFTLD_Astrocytes","C9ALSnoFTLD_Astrocytes","sALSnoFTLD_Astrocytes",
  "control_Microglia","C9noALSnoFTLD_Microglia","C9ALSFTLD_Microglia","C9ALSnoFTLD_Microglia","sALSnoFTLD_Microglia",
  "control_Excitatory","C9noALSnoFTLD_Excitatory","C9ALSFTLD_Excitatory","C9ALSnoFTLD_Excitatory","sALSnoFTLD_Excitatory",
  "control_Inhibitory","C9noALSnoFTLD_Inhibitory","C9ALSFTLD_Inhibitory","C9ALSnoFTLD_Inhibitory","sALSnoFTLD_Inhibitory"
  ))


intATAC$diagnosis_cellsubtype <- paste0(intATAC$diagnosis, "_", intATAC$cellsubtype)

intATAC$diagnosis_cellsubtype<-factor(intATAC$diagnosis_cellsubtype,c(
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


Idents(intATAC)<-"diagnosis"
Idents(intATAC)<-"pseudobulk"
pseudobulk<-c("pseudobulk","pseudobulk","pseudobulk","pseudobulk","pseudobulk")
names(pseudobulk)<-levels(intATAC)
intATAC<-RenameIdents(intATAC,pseudobulk)
intATAC$pseudobulk<-Idents(intATAC)

pseudobulk<-c("snATAC")
names(pseudobulk)<-levels(intATAC)
intATAC<-RenameIdents(intRNA,pseudobulk)
intATAC$pseudobulk<-Idents(intATAC)

saveRDS(intATAC,file="objects/intATAC_final.RDS")

Idents(intATAC)<-"celltype"
intATAC.oligo<-subset(intATAC,idents="Oligodendrocytes")
saveRDS(intATAC.oligo, file="objects/ATAC_oligo.RDS")
intATAC.OPC<-subset(intATAC,idents="OPC")
saveRDS(intATAC.oligo, file="objects/ATAC_OPC.RDS")
intATAC.astro<-subset(intATAC,idents="Astrocytes")
saveRDS(intATAC.oligo, file="objects/ATAC_astrocytes.RDS")
intATAC.micro<-subset(intATAC,idents="Microglia")
saveRDS(intATAC.oligo, file="objects/ATAC_microglia.RDS")
intATAC.ex<-subset(intATAC,idents="excitatory")
saveRDS(intATAC.oligo, file="objects/ATAC_excitatory.RDS")
intATAC.in<-subset(intATAC,idents="inhibitory")
saveRDS(intATAC.oligo, file="objects/ATAC_inhibitory.RDS")

Idents(intATAC)<-"cellsubtype"
intATAC.astro1<-subset(intATAC, idents=c("AST-PP"))
intATAC.astro2<-subset(intATAC, idents=c("AST-FB"))
saveRDS(intATAC.astro1, file="objects/ATAC_AST-PP.RDS")
saveRDS(intATAC.astro2, file="objects/ATAC_AST-FB.RDS")
intATAC.ex1<-subset(intATAC, idents=c("L23"))
intATAC.ex2<-subset(intATAC, idents=c("L4"))
intATAC.ex3<-subset(intATAC, idents=c("L56"))
intATAC.ex4<-subset(intATAC, idents=c(L56-CC))
saveRDS(intATAC.ex1, file="objects/ATAC_L2-3.RDS")
saveRDS(intATAC.ex2, file="objects/ATAC_L4.RDS")
saveRDS(intATAC.ex3, file="objects/ATAC_L5-6.RDS")
saveRDS(intATAC.ex4, file="objects/ATAC_L5-6CC.RDS")
intATAC.in1<-subset(intATAC, idents=c("IN-PV"))
intATAC.in2<-subset(intATAC, idents=c("IN-SST"))
intATAC.in3<-subset(intATAC, idents=c("IN-SV2C"))
intATAC.in4<-subset(intATAC, idents=c("IN-VIP"))
saveRDS(intATAC.in1, file="objects/ATAC_IN-PV.RDS")
saveRDS(intATAC.in2, file="objects/ATAC_IN-SST.RDS")
saveRDS(intATAC.in3, file="objects/ATAC_IN-SV2C.RDS")
saveRDS(intATAC.in4, file="objects/ATAC_IN-VIP.RDS")
