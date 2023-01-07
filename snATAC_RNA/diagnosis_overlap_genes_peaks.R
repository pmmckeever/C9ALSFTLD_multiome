#Celltype Figure analysis (3 onwards)

Idents(intOLIGO)<-'diagnosis'
Idents(intOPC)<-'diagnosis'
Idents(intASTPP)<-'diagnosis'
Idents(intMICRO)<-'diagnosis'

intRNA.OLIGO.avg<-AverageExpression(intOLIGO,assays="RNA",return.seurat=T)
intRNA.OPC.avg<-AverageExpression(intOPC,assays="RNA",return.seurat=T)
intRNA.ASTPP.avg<-AverageExpression(intASTPP,assays="RNA",return.seurat=T)
intRNA.micro.avg<-AverageExpression(intMICRO,assays="RNA",return.seurat=T)

intRNA.OLIGO.avg$diagnosis<-Idents(intRNA.OLIGO.avg)
intRNA.OPC.avg$diagnosis<-Idents(intRNA.OPC.avg)
intRNA.ASTPP.avg$diagnosis<-Idents(intRNA.ASTPP.avg)
intRNA.micro.avg$diagnosis<-Idents(intRNA.micro.avg)


DE_OLIGO<-cellsubtype.C9ALSnoFTLD.cutoff$gene[cellsubtype.C9ALSnoFTLD.cutoff$cellsubtype == "Oligodendrocytes"]
DE_OPC<-cellsubtype.C9ALSnoFTLD.cutoff$gene[cellsubtype.C9ALSnoFTLD.cutoff$cellsubtype == "OPC"]
DE_ASTPP<-cellsubtype.C9ALSnoFTLD.cutoff$gene[cellsubtype.C9ALSnoFTLD.cutoff$cellsubtype == "AST-PP"]
DE_micro<-cellsubtype.C9ALSnoFTLD.cutoff$gene[cellsubtype.C9ALSnoFTLD.cutoff$cellsubtype == "Microglia"]

OLIGO_C9_heatmap<-dittoHeatmap(intRNA.OLIGO.avg,genes=DE_OLIGO[1:50],main = "Oligo DE in C9ALSnoFTLD",
                               annot.by="diagnosis",annot.colors = group_colors, show_colnames = F, show_rownames = T,border_color = "lightgray")
#ggsave('figures/micro_heatmap_C9vsCTRL.pdf', plot=micro_C9_heatmap, width = 5, height = 6, units='in')

OPC_C9_heatmap<-dittoHeatmap(intRNA.OPC.avg,genes=DE_OPC[1:50],main = "OPC DE in C9ALSnoFTLD",
                             annot.by="diagnosis",annot.colors = group_colors, show_colnames = F, show_rownames = T,border_color = "lightgray")
#ggsave('figures/micro_heatmap_C9vsCTRL.pdf', plot=micro_C9_heatmap, width = 5, height = 6, units='in')

ASTPP_C9_heatmap<-dittoHeatmap(intRNA.micro.avg,genes=DE_ASTPP[1:50],main = "ASTPP DE in C9ALSnoFTLD",
                               annot.by="diagnosis",annot.colors = group_colors, show_colnames = F, show_rownames = T,border_color = "lightgray")
#ggsave('figures/micro_heatmap_C9vsCTRL.pdf', plot=micro_C9_heatmap, width = 5, height = 6, units='in')

micro_C9_heatmap<-dittoHeatmap(intRNA.micro.avg,genes=DE_micro[1:50],main = "micro DE in C9ALSnoFTLD",
                               annot.by="diagnosis",annot.colors = group_colors, show_colnames = F, show_rownames = T,border_color = "lightgray")
#ggsave('figures/micro_heatmap_C9vsCTRL.pdf', plot=micro_C9_heatmap, width = 5, height = 6, units='in')

#################################################
###Overlap between peaks and genes C9ALSFTLD
#################################################

load("data/cellsubtype_C9ALSFTLD_markers_LR_RNA_V3only_noRP.rda")
cellsubtype.C9ALSFTLD.gene.markers<-cellsubtype.C9ALSFTLD.markers
cellsubtype.C9ALSFTLD.gene.cutoff<-cellsubtype.C9ALSFTLD.gene.markers[cellsubtype.C9ALSFTLD.gene.markers$p_val_adj<0.01, ]
rm(cellsubtype.C9ALSFTLD.markers)

#load("/mnt/WORKHORSE/Jan22_RNA_ATAC/data/ATAC_differential_peaks/cellsubtype_C9ALSFTLD_diffpeaks_LR.rda")
load("data/cellsubtype_C9ALSFTLD_diffpeaks_LR.rda")
cellsubtype.C9ALSFTLD.peaks.markers<-cellsubtype.C9ALSFTLD.markers
cellsubtype.C9ALSFTLD.peaks.cutoff<-cellsubtype.C9ALSFTLD.markers[cellsubtype.C9ALSFTLD.markers$p_val_adj<0.01, ]
rm(cellsubtype.C9ALSFTLD.markers)

cellsubtype.C9ALSFTLD.gene.cutoff$cellsubtype<-cellsubtype.C9ALSFTLD.gene.cutoff$cellsubtype
cellsubtype.C9ALSFTLD.gene.cutoff$cellsubtype<-gsub("/","",cellsubtype.C9ALSFTLD.gene.cutoff$cellsubtype)
cellsubtype.C9ALSFTLD.peaks.cutoff$cellsubtype<-cellsubtype.C9ALSFTLD.peaks.cutoff$celltype
cellsubtype.C9ALSFTLD.peaks.cutoff$cellsubtype<-gsub("/","",cellsubtype.C9ALSFTLD.peaks.cutoff$cellsubtype)

idents<-unique(cellsubtype.C9ALSFTLD.gene.cutoff$cellsubtype)
overlap_df<-data.frame()
for(i in 1:length(idents)){
  
  cur_cellsubtype <- idents[[i]]
  print(cur_cellsubtype)
  
  genes <- cellsubtype.C9ALSFTLD.gene.cutoff$gene[cellsubtype.C9ALSFTLD.gene.cutoff$cellsubtype == cur_cellsubtype]
  peaks <- cellsubtype.C9ALSFTLD.peaks.cutoff$genes[cellsubtype.C9ALSFTLD.peaks.cutoff$cellsubtype == cur_cellsubtype]
  overlap <- intersect(genes, peaks)
  overlap$cellsubtype<-cur_cellsubtype
  overlap_df <- rbind(overlap_df, overlap)
}
cellsubtype.C9ALSFTLD.markers <- combined_df

#oligo
DE_oligo_genes<-cellsubtype.C9ALSFTLD.gene.cutoff$gene[cellsubtype.C9ALSFTLD.gene.cutoff$cellsubtype == "Oligodendrocytes"]
DE_oligo_peaks<-cellsubtype.C9ALSFTLD.peaks.cutoff$genes[cellsubtype.C9ALSFTLD.peaks.cutoff$cellsubtype == "Oligodendrocytes"]
DE_oligo_overlap<-intersect(DE_oligo_genes,DE_oligo_peaks)
DE_oligo_overlap

#OPC
DE_OPC_genes<-cellsubtype.C9ALSFTLD.gene.cutoff$gene[cellsubtype.C9ALSFTLD.gene.cutoff$cellsubtype == "OPC"]
DE_OPC_peaks<-cellsubtype.C9ALSFTLD.peaks.cutoff$genes[cellsubtype.C9ALSFTLD.peaks.cutoff$cellsubtype == "OPC"]
DE_OPC_overlap<-intersect(DE_OPC_genes,DE_OPC_peaks)
DE_OPC_overlap

#AST-PP
DE_ASTPP_genes<-cellsubtype.C9ALSFTLD.gene.cutoff$gene[cellsubtype.C9ALSFTLD.gene.cutoff$cellsubtype == "AST-PP"]
DE_ASTPP_peaks<-cellsubtype.C9ALSFTLD.peaks.cutoff$genes[cellsubtype.C9ALSFTLD.peaks.cutoff$cellsubtype == "AST-PP"]
DE_ASTPP_overlap<-intersect(DE_ASTPP_genes,DE_ASTPP_peaks)
DE_ASTPP_overlap

#AST-FB
DE_ASTFB_genes<-cellsubtype.C9ALSFTLD.gene.cutoff$gene[cellsubtype.C9ALSFTLD.gene.cutoff$cellsubtype == "AST-FB"]
DE_ASTFB_peaks<-cellsubtype.C9ALSFTLD.peaks.cutoff$genes[cellsubtype.C9ALSFTLD.peaks.cutoff$cellsubtype == "AST-FB"]
DE_ASTFB_overlap<-intersect(DE_ASTFB_genes,DE_ASTFB_peaks)
DE_ASTFB_overlap

#micro
DE_micro_genes<-cellsubtype.C9ALSFTLD.gene.cutoff$gene[cellsubtype.C9ALSFTLD.gene.cutoff$cellsubtype == "Microglia"]
DE_micro_peaks<-cellsubtype.C9ALSFTLD.peaks.cutoff$genes[cellsubtype.C9ALSFTLD.peaks.cutoff$cellsubtype == "Microglia"]
DE_micro_overlap<-intersect(DE_micro_genes,DE_micro_peaks)
DE_micro_overlap

#L23
DE_L23_genes<-cellsubtype.C9ALSFTLD.gene.cutoff$gene[cellsubtype.C9ALSFTLD.gene.cutoff$cellsubtype == "L23"]
DE_L23_peaks<-cellsubtype.C9ALSFTLD.peaks.cutoff$genes[cellsubtype.C9ALSFTLD.peaks.cutoff$cellsubtype == "L23"]
DE_L23_overlap<-intersect(DE_L23_genes,DE_L23_peaks)
DE_L23_overlap

#L4
DE_L4_genes<-cellsubtype.C9ALSFTLD.gene.cutoff$gene[cellsubtype.C9ALSFTLD.gene.cutoff$cellsubtype == "L4"]
DE_L4_peaks<-cellsubtype.C9ALSFTLD.peaks.cutoff$genes[cellsubtype.C9ALSFTLD.peaks.cutoff$cellsubtype == "L4"]
DE_L4_overlap<-intersect(DE_L4_genes,DE_L4_peaks)
DE_L4_overlap

#L56
DE_L56_genes<-cellsubtype.C9ALSFTLD.gene.cutoff$gene[cellsubtype.C9ALSFTLD.gene.cutoff$cellsubtype == "L56"]
DE_L56_peaks<-cellsubtype.C9ALSFTLD.peaks.cutoff$genes[cellsubtype.C9ALSFTLD.peaks.cutoff$cellsubtype == "L56"]
DE_L56_overlap<-intersect(DE_L56_genes,DE_L56_peaks)
DE_L56_overlap

#L56-CC
DE_L56cc_genes<-cellsubtype.C9ALSFTLD.gene.cutoff$gene[cellsubtype.C9ALSFTLD.gene.cutoff$cellsubtype == "L56-CC"]
DE_L56cc_peaks<-cellsubtype.C9ALSFTLD.peaks.cutoff$genes[cellsubtype.C9ALSFTLD.peaks.cutoff$cellsubtype == "L56-CC"]
DE_L56cc_overlap<-intersect(DE_L56cc_genes,DE_L56cc_peaks)
DE_L56cc_overlap

#INPV
DE_INPV_genes<-cellsubtype.C9ALSFTLD.gene.cutoff$gene[cellsubtype.C9ALSFTLD.gene.cutoff$cellsubtype == "IN-PV"]
DE_INPV_peaks<-cellsubtype.C9ALSFTLD.peaks.cutoff$genes[cellsubtype.C9ALSFTLD.peaks.cutoff$cellsubtype == "IN-PV"]
DE_INPV_overlap<-intersect(DE_INPV_genes,DE_INPV_peaks)
DE_INPV_overlap

#INSST
DE_INSST_genes<-cellsubtype.C9ALSFTLD.gene.cutoff$gene[cellsubtype.C9ALSFTLD.gene.cutoff$cellsubtype == "IN-SST"]
DE_INSST_peaks<-cellsubtype.C9ALSFTLD.peaks.cutoff$genes[cellsubtype.C9ALSFTLD.peaks.cutoff$cellsubtype == "IN-SST"]
DE_INSST_overlap<-intersect(DE_INSST_genes,DE_INSST_peaks)
DE_INSST_overlap

#INSV2C
DE_INSV2C_genes<-cellsubtype.C9ALSFTLD.gene.cutoff$gene[cellsubtype.C9ALSFTLD.gene.cutoff$cellsubtype == "IN-SV2C"]
DE_INSV2C_peaks<-cellsubtype.C9ALSFTLD.peaks.cutoff$genes[cellsubtype.C9ALSFTLD.peaks.cutoff$cellsubtype == "IN-SV2C"]
DE_INSV2C_overlap<-intersect(DE_INSV2C_genes,DE_INSV2C_peaks)
DE_INSV2C_overlap

#INVIP
DE_INVIP_genes<-cellsubtype.C9ALSFTLD.gene.cutoff$gene[cellsubtype.C9ALSFTLD.gene.cutoff$cellsubtype == "IN-VIP"]
DE_INVIP_peaks<-cellsubtype.C9ALSFTLD.peaks.cutoff$genes[cellsubtype.C9ALSFTLD.peaks.cutoff$cellsubtype == "IN-VIP"]
DE_INVIP_overlap<-intersect(DE_INVIP_genes,DE_INVIP_peaks)
DE_INVIP_overlap

my_list<-list(DE_oligo_overlap,DE_OPC_overlap,DE_ASTPP_overlap,DE_micro_overlap,DE_L23_overlap,DE_L4_overlap,DE_L56_overlap,DE_L56cc_overlap,DE_INPV_overlap,DE_INSST_overlap,DE_INVIP_overlap)
my_list
names(my_list)<-c("Oligodendrocytes","OPC","AST-PP","Microglia","L23","L4","L56","L56-CC","IN-PV","IN-SST","IN-VIP")
my_list
saveRDS(my_list,file="overlap_DEG_DAR_C9ALSFTLD.rds")
write.xlsx(my_list,file="overlap_DEG_DAR_C9ALSFTLD.xlsx")

#################################################
###Overlap between peaks and genes C9ALSnoFTLD
#################################################

load("/mnt/WORKHORSE/Jan22_RNA_ATAC/V3only/cellsubtype_C9ALSnoFTLD_markers_LR_RNA_V3only_noRP.rda")
cellsubtype.C9ALSnoFTLD.gene.markers<-cellsubtype.C9ALSnoFTLD.markers
cellsubtype.C9ALSnoFTLD.gene.cutoff<-cellsubtype.C9ALSnoFTLD.gene.markers[cellsubtype.C9ALSnoFTLD.gene.markers$p_val_adj<0.01, ]
rm(cellsubtype.C9ALSnoFTLD.markers)

#load("/mnt/WORKHORSE/Jan22_RNA_ATAC/data/ATAC_differential_peaks/cellsubtype_C9ALSnoFTLD_diffpeaks_LR.rda")
load("/mnt/WORKHORSE/Jan22_RNA_ATAC/tables/closest_genes/cellsubtype_C9ALSnoFTLD_diffpeaks_LR.rda")
cellsubtype.C9ALSnoFTLD.peaks.markers<-cellsubtype.C9ALSnoFTLD.markers
cellsubtype.C9ALSnoFTLD.peaks.cutoff<-cellsubtype.C9ALSnoFTLD.markers[cellsubtype.C9ALSnoFTLD.markers$p_val_adj<0.01, ]
rm(cellsubtype.C9ALSnoFTLD.markers)

cellsubtype.C9ALSnoFTLD.gene.cutoff$cellsubtype<-cellsubtype.C9ALSnoFTLD.gene.cutoff$cellsubtype
cellsubtype.C9ALSnoFTLD.gene.cutoff$cellsubtype<-gsub("/","",cellsubtype.C9ALSnoFTLD.gene.cutoff$cellsubtype)
cellsubtype.C9ALSnoFTLD.peaks.cutoff$cellsubtype<-cellsubtype.C9ALSnoFTLD.peaks.cutoff$celltype
cellsubtype.C9ALSnoFTLD.peaks.cutoff$cellsubtype<-gsub("/","",cellsubtype.C9ALSnoFTLD.peaks.cutoff$cellsubtype)

#oligo
DE_oligo_genes<-cellsubtype.C9ALSnoFTLD.gene.cutoff$gene[cellsubtype.C9ALSnoFTLD.gene.cutoff$cellsubtype == "Oligodendrocytes"]
DE_oligo_peaks<-cellsubtype.C9ALSnoFTLD.peaks.cutoff$genes[cellsubtype.C9ALSnoFTLD.peaks.cutoff$cellsubtype == "Oligodendrocytes"]
DE_oligo_overlap<-intersect(DE_oligo_genes,DE_oligo_peaks)
DE_oligo_overlap

#OPC
DE_OPC_genes<-cellsubtype.C9ALSnoFTLD.gene.cutoff$gene[cellsubtype.C9ALSnoFTLD.gene.cutoff$cellsubtype == "OPC"]
DE_OPC_peaks<-cellsubtype.C9ALSnoFTLD.peaks.cutoff$genes[cellsubtype.C9ALSnoFTLD.peaks.cutoff$cellsubtype == "OPC"]
DE_OPC_overlap<-intersect(DE_OPC_genes,DE_OPC_peaks)
DE_OPC_overlap

#AST-PP
DE_ASTPP_genes<-cellsubtype.C9ALSnoFTLD.gene.cutoff$gene[cellsubtype.C9ALSnoFTLD.gene.cutoff$cellsubtype == "AST-PP"]
DE_ASTPP_peaks<-cellsubtype.C9ALSnoFTLD.peaks.cutoff$genes[cellsubtype.C9ALSnoFTLD.peaks.cutoff$cellsubtype == "AST-PP"]
DE_ASTPP_overlap<-intersect(DE_ASTPP_genes,DE_ASTPP_peaks)
DE_ASTPP_overlap

#AST-FB
DE_ASTFB_genes<-cellsubtype.C9ALSnoFTLD.gene.cutoff$gene[cellsubtype.C9ALSnoFTLD.gene.cutoff$cellsubtype == "AST-FB"]
DE_ASTFB_peaks<-cellsubtype.C9ALSnoFTLD.peaks.cutoff$genes[cellsubtype.C9ALSnoFTLD.peaks.cutoff$cellsubtype == "AST-FB"]
DE_ASTFB_overlap<-intersect(DE_ASTFB_genes,DE_ASTFB_peaks)
DE_ASTFB_overlap

#micro
DE_micro_genes<-cellsubtype.C9ALSnoFTLD.gene.cutoff$gene[cellsubtype.C9ALSnoFTLD.gene.cutoff$cellsubtype == "Microglia"]
DE_micro_peaks<-cellsubtype.C9ALSnoFTLD.peaks.cutoff$genes[cellsubtype.C9ALSnoFTLD.peaks.cutoff$cellsubtype == "Microglia"]
DE_micro_overlap<-intersect(DE_micro_genes,DE_micro_peaks)
DE_micro_overlap

#L23
DE_L23_genes<-cellsubtype.C9ALSnoFTLD.gene.cutoff$gene[cellsubtype.C9ALSnoFTLD.gene.cutoff$cellsubtype == "L23"]
DE_L23_peaks<-cellsubtype.C9ALSnoFTLD.peaks.cutoff$genes[cellsubtype.C9ALSnoFTLD.peaks.cutoff$cellsubtype == "L23"]
DE_L23_overlap<-intersect(DE_L23_genes,DE_L23_peaks)
DE_L23_overlap

#L4
DE_L4_genes<-cellsubtype.C9ALSnoFTLD.gene.cutoff$gene[cellsubtype.C9ALSnoFTLD.gene.cutoff$cellsubtype == "L4"]
DE_L4_peaks<-cellsubtype.C9ALSnoFTLD.peaks.cutoff$genes[cellsubtype.C9ALSnoFTLD.peaks.cutoff$cellsubtype == "L4"]
DE_L4_overlap<-intersect(DE_L4_genes,DE_L4_peaks)
DE_L4_overlap

#L56
DE_L56_genes<-cellsubtype.C9ALSnoFTLD.gene.cutoff$gene[cellsubtype.C9ALSnoFTLD.gene.cutoff$cellsubtype == "L56"]
DE_L56_peaks<-cellsubtype.C9ALSnoFTLD.peaks.cutoff$genes[cellsubtype.C9ALSnoFTLD.peaks.cutoff$cellsubtype == "L56"]
DE_L56_overlap<-intersect(DE_L56_genes,DE_L56_peaks)
DE_L56_overlap

#L56-CC
DE_L56cc_genes<-cellsubtype.C9ALSnoFTLD.gene.cutoff$gene[cellsubtype.C9ALSnoFTLD.gene.cutoff$cellsubtype == "L56-CC"]
DE_L56cc_peaks<-cellsubtype.C9ALSnoFTLD.peaks.cutoff$genes[cellsubtype.C9ALSnoFTLD.peaks.cutoff$cellsubtype == "L56-CC"]
DE_L56cc_overlap<-intersect(DE_L56cc_genes,DE_L56cc_peaks)
DE_L56cc_overlap

#INPV
DE_INPV_genes<-cellsubtype.C9ALSnoFTLD.gene.cutoff$gene[cellsubtype.C9ALSnoFTLD.gene.cutoff$cellsubtype == "IN-PV"]
DE_INPV_peaks<-cellsubtype.C9ALSnoFTLD.peaks.cutoff$genes[cellsubtype.C9ALSnoFTLD.peaks.cutoff$cellsubtype == "IN-PV"]
DE_INPV_overlap<-intersect(DE_INPV_genes,DE_INPV_peaks)
DE_INPV_overlap

#INSST
DE_INSST_genes<-cellsubtype.C9ALSnoFTLD.gene.cutoff$gene[cellsubtype.C9ALSnoFTLD.gene.cutoff$cellsubtype == "IN-SST"]
DE_INSST_peaks<-cellsubtype.C9ALSnoFTLD.peaks.cutoff$genes[cellsubtype.C9ALSnoFTLD.peaks.cutoff$cellsubtype == "IN-SST"]
DE_INSST_overlap<-intersect(DE_INSST_genes,DE_INSST_peaks)
DE_INSST_overlap

#INSV2C
DE_INSV2C_genes<-cellsubtype.C9ALSnoFTLD.gene.cutoff$gene[cellsubtype.C9ALSnoFTLD.gene.cutoff$cellsubtype == "IN-SV2C"]
DE_INSV2C_peaks<-cellsubtype.C9ALSnoFTLD.peaks.cutoff$genes[cellsubtype.C9ALSnoFTLD.peaks.cutoff$cellsubtype == "IN-SV2C"]
DE_INSV2C_overlap<-intersect(DE_INSV2C_genes,DE_INSV2C_peaks)
DE_INSV2C_overlap

#INVIP
DE_INVIP_genes<-cellsubtype.C9ALSnoFTLD.gene.cutoff$gene[cellsubtype.C9ALSnoFTLD.gene.cutoff$cellsubtype == "IN-VIP"]
DE_INVIP_peaks<-cellsubtype.C9ALSnoFTLD.peaks.cutoff$genes[cellsubtype.C9ALSnoFTLD.peaks.cutoff$cellsubtype == "IN-VIP"]
DE_INVIP_overlap<-intersect(DE_INVIP_genes,DE_INVIP_peaks)
DE_INVIP_overlap


my_list<-list(DE_oligo_overlap,DE_OPC_overlap,DE_ASTPP_overlap,DE_micro_overlap,DE_L23_overlap,DE_L4_overlap,DE_L56cc_overlap,DE_INPV_overlap,DE_INSST_overlap,DE_INVIP_overlap)
my_list
names(my_list)<-c("Oligodendrocytes","OPC","AST-PP","Microglia","L23","L4","L56-CC","IN-PV","IN-SST","IN-VIP")
my_list
saveRDS(my_list,file="overlap_DEG_DAR_C9ALSnoFTLD.rds")
write.xlsx(my_list,file="overlap_DEG_DAR_C9ALSnoFTLD.xlsx")

#################################################
###Overlap between peaks and genes sALSnoFTLD
#################################################

load("/mnt/WORKHORSE/Jan22_RNA_ATAC/V3only/cellsubtype_sALSnoFTLD_markers_LR_RNA_V3only_noRP.rda")
cellsubtype.sALSnoFTLD.gene.markers<-cellsubtype.sALSnoFTLD.markers
cellsubtype.sALSnoFTLD.gene.cutoff<-cellsubtype.sALSnoFTLD.gene.markers[cellsubtype.sALSnoFTLD.gene.markers$p_val_adj<0.01, ]
rm(cellsubtype.sALSnoFTLD.markers)

#load("/mnt/WORKHORSE/Jan22_RNA_ATAC/data/ATAC_differential_peaks/cellsubtype_sALSnoFTLD_diffpeaks_LR.rda")
load("/mnt/WORKHORSE/Jan22_RNA_ATAC/tables/closest_genes/cellsubtype_sALSnoFTLD_diffpeaks_LR.rda")
cellsubtype.sALSnoFTLD.peaks.markers<-cellsubtype.sALSnoFTLD.markers
cellsubtype.sALSnoFTLD.peaks.cutoff<-cellsubtype.sALSnoFTLD.markers[cellsubtype.sALSnoFTLD.markers$p_val_adj<0.01, ]
rm(cellsubtype.sALSnoFTLD.markers)

cellsubtype.sALSnoFTLD.gene.cutoff$cellsubtype<-cellsubtype.sALSnoFTLD.gene.cutoff$cellsubtype
cellsubtype.sALSnoFTLD.gene.cutoff$cellsubtype<-gsub("/","",cellsubtype.sALSnoFTLD.gene.cutoff$cellsubtype)
cellsubtype.sALSnoFTLD.peaks.cutoff$cellsubtype<-cellsubtype.sALSnoFTLD.peaks.cutoff$celltype
cellsubtype.sALSnoFTLD.peaks.cutoff$cellsubtype<-gsub("/","",cellsubtype.sALSnoFTLD.peaks.cutoff$cellsubtype)

#oligo
DE_oligo_genes<-cellsubtype.sALSnoFTLD.gene.cutoff$gene[cellsubtype.sALSnoFTLD.gene.cutoff$cellsubtype == "Oligodendrocytes"]
DE_oligo_peaks<-cellsubtype.sALSnoFTLD.peaks.cutoff$genes[cellsubtype.sALSnoFTLD.peaks.cutoff$cellsubtype == "Oligodendrocytes"]
DE_oligo_overlap<-intersect(DE_oligo_genes,DE_oligo_peaks)
DE_oligo_overlap

#OPC
DE_OPC_genes<-cellsubtype.sALSnoFTLD.gene.cutoff$gene[cellsubtype.sALSnoFTLD.gene.cutoff$cellsubtype == "OPC"]
DE_OPC_peaks<-cellsubtype.sALSnoFTLD.peaks.cutoff$genes[cellsubtype.sALSnoFTLD.peaks.cutoff$cellsubtype == "OPC"]
DE_OPC_overlap<-intersect(DE_OPC_genes,DE_OPC_peaks)
DE_OPC_overlap

#AST-PP
DE_ASTPP_genes<-cellsubtype.sALSnoFTLD.gene.cutoff$gene[cellsubtype.sALSnoFTLD.gene.cutoff$cellsubtype == "AST-PP"]
DE_ASTPP_peaks<-cellsubtype.sALSnoFTLD.peaks.cutoff$genes[cellsubtype.sALSnoFTLD.peaks.cutoff$cellsubtype == "AST-PP"]
DE_ASTPP_overlap<-intersect(DE_ASTPP_genes,DE_ASTPP_peaks)
DE_ASTPP_overlap

#AST-FB
DE_ASTFB_genes<-cellsubtype.sALSnoFTLD.gene.cutoff$gene[cellsubtype.sALSnoFTLD.gene.cutoff$cellsubtype == "AST-FB"]
DE_ASTFB_peaks<-cellsubtype.sALSnoFTLD.peaks.cutoff$genes[cellsubtype.sALSnoFTLD.peaks.cutoff$cellsubtype == "AST-FB"]
DE_ASTFB_overlap<-intersect(DE_ASTFB_genes,DE_ASTFB_peaks)
DE_ASTFB_overlap

#micro
DE_micro_genes<-cellsubtype.sALSnoFTLD.gene.cutoff$gene[cellsubtype.sALSnoFTLD.gene.cutoff$cellsubtype == "Microglia"]
DE_micro_peaks<-cellsubtype.sALSnoFTLD.peaks.cutoff$genes[cellsubtype.sALSnoFTLD.peaks.cutoff$cellsubtype == "Microglia"]
DE_micro_overlap<-intersect(DE_micro_genes,DE_micro_peaks)
DE_micro_overlap

#L23
DE_L23_genes<-cellsubtype.sALSnoFTLD.gene.cutoff$gene[cellsubtype.sALSnoFTLD.gene.cutoff$cellsubtype == "L23"]
DE_L23_peaks<-cellsubtype.sALSnoFTLD.peaks.cutoff$genes[cellsubtype.sALSnoFTLD.peaks.cutoff$cellsubtype == "L23"]
DE_L23_overlap<-intersect(DE_L23_genes,DE_L23_peaks)
DE_L23_overlap

#L4
DE_L4_genes<-cellsubtype.sALSnoFTLD.gene.cutoff$gene[cellsubtype.sALSnoFTLD.gene.cutoff$cellsubtype == "L4"]
DE_L4_peaks<-cellsubtype.sALSnoFTLD.peaks.cutoff$genes[cellsubtype.sALSnoFTLD.peaks.cutoff$cellsubtype == "L4"]
DE_L4_overlap<-intersect(DE_L4_genes,DE_L4_peaks)
DE_L4_overlap

#L56
DE_L56_genes<-cellsubtype.sALSnoFTLD.gene.cutoff$gene[cellsubtype.sALSnoFTLD.gene.cutoff$cellsubtype == "L56"]
DE_L56_peaks<-cellsubtype.sALSnoFTLD.peaks.cutoff$genes[cellsubtype.sALSnoFTLD.peaks.cutoff$cellsubtype == "L56"]
DE_L56_overlap<-intersect(DE_L56_genes,DE_L56_peaks)
DE_L56_overlap

#L56-CC
DE_L56cc_genes<-cellsubtype.sALSnoFTLD.gene.cutoff$gene[cellsubtype.sALSnoFTLD.gene.cutoff$cellsubtype == "L56-CC"]
DE_L56cc_peaks<-cellsubtype.sALSnoFTLD.peaks.cutoff$genes[cellsubtype.sALSnoFTLD.peaks.cutoff$cellsubtype == "L56-CC"]
DE_L56cc_overlap<-intersect(DE_L56cc_genes,DE_L56cc_peaks)
DE_L56cc_overlap

#INPV
DE_INPV_genes<-cellsubtype.sALSnoFTLD.gene.cutoff$gene[cellsubtype.sALSnoFTLD.gene.cutoff$cellsubtype == "IN-PV"]
DE_INPV_peaks<-cellsubtype.sALSnoFTLD.peaks.cutoff$genes[cellsubtype.sALSnoFTLD.peaks.cutoff$cellsubtype == "IN-PV"]
DE_INPV_overlap<-intersect(DE_INPV_genes,DE_INPV_peaks)
DE_INPV_overlap

#INSST
DE_INSST_genes<-cellsubtype.sALSnoFTLD.gene.cutoff$gene[cellsubtype.sALSnoFTLD.gene.cutoff$cellsubtype == "IN-SST"]
DE_INSST_peaks<-cellsubtype.sALSnoFTLD.peaks.cutoff$genes[cellsubtype.sALSnoFTLD.peaks.cutoff$cellsubtype == "IN-SST"]
DE_INSST_overlap<-intersect(DE_INSST_genes,DE_INSST_peaks)
DE_INSST_overlap

#INSV2C
DE_INSV2C_genes<-cellsubtype.sALSnoFTLD.gene.cutoff$gene[cellsubtype.sALSnoFTLD.gene.cutoff$cellsubtype == "IN-SV2C"]
DE_INSV2C_peaks<-cellsubtype.sALSnoFTLD.peaks.cutoff$genes[cellsubtype.sALSnoFTLD.peaks.cutoff$cellsubtype == "IN-SV2C"]
DE_INSV2C_overlap<-intersect(DE_INSV2C_genes,DE_INSV2C_peaks)
DE_INSV2C_overlap

#INVIP
DE_INVIP_genes<-cellsubtype.sALSnoFTLD.gene.cutoff$gene[cellsubtype.sALSnoFTLD.gene.cutoff$cellsubtype == "IN-VIP"]
DE_INVIP_peaks<-cellsubtype.sALSnoFTLD.peaks.cutoff$genes[cellsubtype.sALSnoFTLD.peaks.cutoff$cellsubtype == "IN-VIP"]
DE_INVIP_overlap<-intersect(DE_INVIP_genes,DE_INVIP_peaks)
DE_INVIP_overlap

my_list<-list(DE_oligo_overlap,DE_OPC_overlap,DE_ASTPP_overlap,DE_micro_overlap,DE_L23_overlap,DE_L4_overlap,DE_INPV_overlap,DE_INSST_overlap,DE_INVIP_overlap)
my_list
names(my_list)<-c("Oligodendrocytes","OPC","AST-PP","Microglia","L23","L4","IN-PV","IN-SST","IN-VIP")
my_list
saveRDS(my_list,file="overlap_DEG_DAR_sALSnoFTLD.rds")
write.xlsx(my_list,file="overlap_DEG_DAR_sALSnoFTLD.xlsx")

# load ggvenn packag
elibrary("ggvenn")

# use list as input
DE_micro_venn <-list('DE genes'=c(DE_micro_genes),'DE peaks'=c(DE_micro_peaks))

# creating venn diagram for four sets
# and displaying only two sets
ggvenn(DE_micro_venn,c("DE genes","DE peaks"),show_percentage=FALSE,
       fill_color=c("","orange"))

fragments<-CreateFragmentObject(path="/mnt/WORKHORSE/Jan22_RNA_ATAC/fragments/fragments.tsv.gz",cells=colnames(intL23))
#intATAC<-readRDS(file="~/projects/def-rogaeva/prime235/Oct21_ATAC/objects/ALS_intATAC_final.RDS")
#fragments<-CreateFragmentObject(path="~/projects/def-rogaeva/prime235/Oct21_ATAC/fragments/fragments.tsv.gz",cells=colnames(intATAC))
Fragments(intL23)<-NULL
Fragments(intL23)<-fragments

fragments<-CreateFragmentObject(path="/mnt/WORKHORSE/Jan22_RNA_ATAC/fragments/fragments.tsv.gz",cells=colnames(intL4))
#intATAC<-readRDS(file="~/projects/def-rogaeva/prime235/Oct21_ATAC/objects/ALS_intATAC_final.RDS")
#fragments<-CreateFragmentObject(path="~/projects/def-rogaeva/prime235/Oct21_ATAC/fragments/fragments.tsv.gz",cells=colnames(intATAC))
Fragments(intL4)<-NULL
Fragments(intL4)<-fragments