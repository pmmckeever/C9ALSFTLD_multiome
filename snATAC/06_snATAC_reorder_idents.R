library(Seurat)
library(Signac)
set.seed(1234)

plan("multicore",workers=32)
options(future.globals.maxSize = 125 * 1024^3, seed=TRUE)

intATAC<-readRDS(file="objects/Nov17_labelxfer.RDS")

Idents(intATAC)<-"cellsubtype"
table(intATAC$cellsubtype)

intATAC$cellsubtype <- factor(intATAC$cellsubtype, levels = c("Oligodendrocytes","OPC","AST-PP","AST-FB","Microglia","Endothelial",
                                                                "L2/3","L4","L5/6","L5/6-CC","IN-PV","IN-SST","IN-SV2C","IN-VIP"))

intATAC$diagnosis_cellsubtype <- paste0(intATAC$diagnosis, "_", intATAC$cellsubtype)

Idents(intATAC)<-"cellsubtype"
celltypes<-c("Oligodendrocytes","OPC","Astrocytes","Astrocytes",
    "Microglia","Endothelial","Excitatory","Excitatory","Excitatory",
    "Excitatory","Inhibitory","Inhibitory","Inhibitory","Inhibitory")
names(celltypes)<-levels(intATAC)
intATAC<-RenameIdents(intATAC,celltypes)
intATAC$celltype<-Idents(intATAC)
intATAC$diagnosis_celltype <- paste0(intATAC$diagnosis, "_", intATAC$celltype)

saveRDS(intATAC,file="objects/ATAC_annotated.RDS")