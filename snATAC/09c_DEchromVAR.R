library(Seurat)
library(Signac)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dittoSeq)
library(scCustomize)
library(future)
set.seed(1234)

intEXC <- readRDS("/mnt/WORKHORSE/Jan22_RNA_ATAC/objects_ATAC/Mar11_ATAC_excitatory.RDS")
intOLIGO <- readRDS("/mnt/WORKHORSE/Jan22_RNA_ATAC/objects_ATAC/Mar11_ATAC_oligo.RDS")
intINH <- readRDS("/mnt/WORKHORSE/Jan22_RNA_ATAC/objects_ATAC/Mar11_ATAC_inhibitory.RDS")
intAST <- readRDS("/mnt/WORKHORSE/Jan22_RNA_ATAC/objects_ATAC/Mar11_ATAC_astrocytes.RDS")
intOPC <- readRDS("/mnt/WORKHORSE/Jan22_RNA_ATAC/objects_ATAC/Mar11_ATAC_OPC.RDS")
intMICRO <- readRDS("/mnt/WORKHORSE/Jan22_RNA_ATAC/objects_ATAC/Mar11_ATAC_microglia.RDS")
avg_chromvar_exc <- readRDS("/mnt/WORKHORSE/Jan22_RNA_ATAC/objects_ATAC/Mar11_avg_chromvar_excitatory.RDS")
avg_chromvar_oligo <- readRDS("/mnt/WORKHORSE/Jan22_RNA_ATAC/objects_ATAC/Mar11_avg_chromvar_oligo.RDS")
avg_chromvar_inh <- readRDS("/mnt/WORKHORSE/Jan22_RNA_ATAC/objects_ATAC/Mar11_avg_chromvar_inhibitory.RDS")
avg_chromvar_astro <- readRDS("/mnt/WORKHORSE/Jan22_RNA_ATAC/objects_ATAC/Mar11_avg_chromvar_astro.RDS")
avg_chromvar_micro <- readRDS("/mnt/WORKHORSE/Jan22_RNA_ATAC/objects_ATAC/Mar11_avg_chromvar_micro.RDS")
avg_chromvar_OPC <- readRDS("/mnt/WORKHORSE/Jan22_RNA_ATAC/objects_ATAC/Mar11_avg_chromvar_OPC.RDS")
load("/mnt/WORKHORSE/Jan22_RNA_ATAC/tables/motif_analysis/inh.motif.ctrlvC9ALSFTLD.rda")
load("/mnt/WORKHORSE/Jan22_RNA_ATAC/tables/motif_analysis/astro.motif.ctrlvC9ALSFTLD.rda")
load("/mnt/WORKHORSE/Jan22_RNA_ATAC/tables/motif_analysis/exc.motif.ctrlvC9ALSFTLD.rda")
load("/mnt/WORKHORSE/Jan22_RNA_ATAC/tables/motif_analysis/micro.motif.ctrlvC9ALSFTLD.rda")
load("/mnt/WORKHORSE/Jan22_RNA_ATAC/tables/motif_analysis/oligo.motif.ctrlvC9ALSFTLD.rda")
load("/mnt/WORKHORSE/Jan22_RNA_ATAC/tables/motif_analysis/OPC.motif.ctrlvC9ALSFTLD.rda")

intOLIGO$diagnosis<-factor(intOLIGO$diagnosis,levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))
intOPC$diagnosis<-factor(intOPC$diagnosis,levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))
intAST$diagnosis<-factor(intAST$diagnosis,levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))
intMICRO$diagnosis<-factor(intMICRO$diagnosis,levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))
intEXC$diagnosis<-factor(intEXC$diagnosis,levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))
intINH$diagnosis<-factor(intINH$diagnosis,levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))

avg_chromvar_oligo$diagnosis<-factor(avg_chromvar_oligo$diagnosis, levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))
avg_chromvar_OPC$diagnosis<-factor(avg_chromvar_OPC$diagnosis, levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))
avg_chromvar_astro$diagnosis<-factor(avg_chromvar_astro$diagnosis, levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))
avg_chromvar_micro$diagnosis<-factor(avg_chromvar_micro$diagnosis, levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))
avg_chromvar_exc$diagnosis<-factor(avg_chromvar_exc$diagnosis, levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))
avg_chromvar_inh$diagnosis<-factor(avg_chromvar_inh$diagnosis, levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))


diagnosis_colors<-c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2")
diagnosis<-c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD")

celltype_diagnosis_colors<-c("#E69F00","#009F73","#F0E442","#0072B2","#D55E00","#CC79A7","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2")

,
"#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
"#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
"#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
"#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
"#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
OLIGO OPC ASTO MIC EX IN

#read and add fragment path
intATAC <- readRDS("objects/intATAC_final.RDS")
fragments<-CreateFragmentObject(path="fragments/fragments.tsv.gz",cells=colnames(intATAC))
DefaultAssay(intATAC)<-"peaks"
Fragments(intATAC)<-NULL
Fragments(intATAC)<-fragments


intATAC2<-DietSeurat(intATAC,counts=T,data=T,scale.data=T,assays=c("peaks","RNA","chromvar"),dimreducs=c("lsi", "umap", "harmony"))

#Seurat not transfering harmony and lsi, so move embeddings manually
intATAC2@reductions$harmony <- intATAC@reductions$harmony
intATAC2@reductions$lsi <- intATAC@reductions$lsi

intATAC<-intATAC2
rm(intATAC2)
gc()

DefaultAssay(intATAC)<-'chromvar'
Idents(intATAC)<-'celltype'
celltype_motif_markers<-FindAllMarkers(intATAC)
save(celltype_motif_markers,file="tables/celltype_motif_markers.rda")
write.csv(celltype_motif_markers,file="tables/celltype_motif_markers.csv")

#SET IDENT NAMES
DefaultAssay(intATAC2)<-"chromvar"

intATAC2$diagnosis<-intATAC2$diagnosis
Idents(intATAC2)<-"diagnosis"
intATAC2$diagnosis_celltype <- paste0(intATAC2$celltype, "_", intATAC2$diagnosis)

Idents(intATAC2)<-'diagnosis_celltype'
intATAC2$diagnosis_celltype<-factor(intATAC2$diagnosis_celltype,levels=c(
	"Oligodendrocytes_control","Oligodendrocytes_C9noALSnoFTLD","Oligodendrocytes_C9ALSFTLD","Oligodendrocytes_C9ALSnoFTLD","Oligodendrocytes_sALSnoFTLD",
	"OPC_control","OPC_C9noALSnoFTLD","OPC_C9ALSFTLD","OPC_C9ALSnoFTLD","OPC_sALSnoFTLD",
	"Astrocytes_control","Astrocytes_C9noALSnoFTLD","Astrocytes_C9ALSFTLD","Astrocytes_C9ALSnoFTLD","Astrocytes_sALSnoFTLD",
	"Microglia_control","Microglia_C9noALSnoFTLD","Microglia_C9ALSFTLD","Microglia_C9ALSnoFTLD","Microglia_sALSnoFTLD",
	"Excitatory_control","Excitatory_C9noALSnoFTLD","Excitatory_C9ALSFTLD","Excitatory_C9ALSnoFTLD","Excitatory_sALSnoFTLD",
	"Inhibitory_control","Inhibitory_C9noALSnoFTLD","Inhibitory_C9ALSFTLD","Inhibitory_C9ALSnoFTLD","Inhibitory_sALSnoFTLD"))

avg_chromvar_ATAC<-AverageExpression(intATAC,assays="chromvar",return.seurat=T)
Idents(avg_chromvar_ATAC)
avg_chromvar_ATAC$celltype_diagnosis<-Idents(avg_chromvar_ATAC)
diagnoses<-c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD",
"control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD",
"control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD",
"control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD",
"control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD",
"control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD")
names(diagnoses)<-levels(avg_chromvar_ATAC)
avg_chromvar_ATAC<-RenameIdents(avg_chromvar_ATAC,diagnoses)
avg_chromvar_ATAC$diagnoses<-Idents(avg_chromvar_ATAC)
avg_chromvar_ATAC$diagnoses
avg_chromvar_ATAC$orig.ident
avg_chromvar_ATAC$Celltype<-avg_chromvar_ATAC$orig.ident
avg_chromvar_ATAC$Celltype<-factor(avg_chromvar_ATAC$Celltype,levels=c("Oligodendrocytes","OPC","Astrocytes","Microglia","Excitatory","Inhibitory"))

celltype_motif_heatmap<-dittoHeatmap(avg_chromvar_ATAC,genes=top5,annot.by=c("Celltype","diagnoses"),order.by="Celltype",show_colnames = F, show_rownames = T,cluster_rows=F,complex=F,
annot.colors=celltype_diagnosis_colors, name="z-score",
border_color="lightgray", scaled.to.max=T,heatmap.colors.max.scaled = viridis(50, direction = 1,option="D"),range=seq(-3,3,length.out=51))
ggsave('celltype_motif_heatmap.png',celltype_motif_heatmap)

OPC.motif.ctrlvC9ALSFTLD
oligo.motif.ctrlvC9ALSFTLD
exc.motif.ctrlvC9ALSFTLD
inh.motif.ctrlvC9ALSFTLD
micro.motif.ctrlvC9ALSFTLD
astro.motif.ctrlvC9ALSFTLD


head(oligo.motif.ctrlvC9ALSFTLD)
head(OPC.motif.ctrlvC9ALSFTLD)
head(astro.motif.ctrlvC9ALSFTLD)
head(micro.motif.ctrlvC9ALSFTLD)
head(exc.motif.ctrlvC9ALSFTLD)
head(inh.motif.ctrlvC9ALSFTLD)


celltype_diagnosis_colors<-c("#D55E00","#CC79A7","#E69F00","#0072B2","#009F73","#F0E442","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2")
celltype_colors<-c("#D55E00","#CC79A7","#E69F00","#0072B2","#009F73","#F0E442")
diagnosis_colors<-c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2")

ATAC_chromvar_C9vsCTRL<-dittoHeatmap(avg_chromvar_ATAC,genes=rownames(diff.act.ATAC)[1:100],main = "ATAC DE motifs in C9ALSFTLD",
             annot.by="diagnosis",annot.colors = diagnosis_colors, show_colnames = F, show_rownames = F)
ggsave('figures/ATAC_chromvar_C9vsCTRL.png', plot=ATAC_chromvar_C9vsCTRL, width = 5, height = 6, units='in')


top10<-Extract_Top_Markers(
celltype_motif_markers,
num_genes = 10,
diagnosis_by = "cluster",
rank_by = "avg_log2FC",
gene_column = "gene",
gene_rownames_to_column = FALSE,
data_frame = FALSE,
named_vector = TRUE,
make_unique = FALSE
)

top5<-Extract_Top_Markers(
  celltype_motif_markers,
  num_genes = 5,
  diagnosis_by = "cluster",
  rank_by = "avg_log2FC",
  gene_column = "gene",
  gene_rownames_to_column = FALSE,
  data_frame = FALSE,
  named_vector = TRUE,
  make_unique = FALSE
)


dittoHeatmap(avg_chromvar_ATAC,genes=top5,main = "ATAC DE motifs by celltype",
             annot.by=c("Celltype","diagnoses"),genes=rownames(top5),
             show_colnames = F, show_rownames = T,cluster_rows=T,complex=F,
             annot.colors=diagnosis_colors, name="z-score", 
             border_color="lightgray", heatmap.colors.max.scaled = viridis(50, direction = 1),range=seq(-5,5,length.out=51))

avg_chromvar_astro<-readRDS(file="objects_ATAC/Mar11_avg_chromvar_astro.RDS")
avg_chromvar_micro<-readRDS(file="objects_ATAC/Mar11_avg_chromvar_micro.RDS")
avg_chromvar_excitatory<-readRDS(file="objects_ATAC/Mar11_avg_chromvar_excitatory.RDS")
avg_chromvar_inhibitory<-readRDS(file="objects_ATAC/Mar11_avg_chromvar_inhibitory.RDS")
avg_chromvar_oligo<-readRDS(file="objects_ATAC/Mar11_avg_chromvar_oligo.RDS")
avg_chromvar_OPC<-readRDS(file="objects_ATAC/Mar11_avg_chromvar_OPC.RDS")

oligo.motif.ctrlvC9ALSFTLD
OPC.motif.ctrlvC9ALSFTLD
astro.motif.ctrlvC9ALSFTLD
micro.motif.ctrlvC9ALSFTLD
exc.motif.ctrlvC9ALSFTLD
inh.motif.ctrlvC9ALSFTLD

## For pheatmap_1.0.8 and later:
draw_colnames_45 <- function (coln, gaps, ...) {
    coord = pheatmap:::find_coordinates(length(coln), gaps)
    x = coord$coord - 0.5 * coord$size
    res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
    return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
ns=asNamespace("pheatmap"))

oligo_chromvar_C9vsCTRL<-dittoHeatmap(avg_chromvar_oligo,genes=rownames(oligo.motif.ctrlvC9ALSFTLD)[1:50],
                                      show_colnames = F, show_rownames = T,cluster_rows=T,complex=F,
                                      annot.colors=diagnosis_colors, name="z-score",annot.by="diagnosis", 
                                      border_color="lightgray",treeheight_row=20,fontsize=10,show_legend=F,order.by="diagnosis")

oligo_chromvar_C9vsCTRL<-dittoHeatmap(avg_chromvar_oligo,genes=rownames(oligo.motif.ctrlvC9ALSFTLD)[1:50],
             show_colnames = T, show_rownames = T,cluster_rows=T,complex=F,
             annot.colors=diagnosis_colors, name="z-score", 
             border_color="lightgray",treeheight_row=30,fontsize=7,fontsize_col=8,show_legend=T,order.by="diagnosis")
ggsave('figures/oligo_chromvar_C9vsCTRL.png', plot=oligo_chromvar_C9vsCTRL, width = 2.5, height = 6, units='in')

oligo_motif<-MotifPlot(
  object = intOLIGO,
  motifs = head(rownames(oligo.motif.ctrlvC9ALSFTLD),n=10),
  assay = 'peaks')
ggsave('figures/oligo_motif_C9vsCTRL.png',plot=oligo_motif,width=8, height=6, units='in')

OPC_chromvar_C9vsCTRL<-dittoHeatmap(avg_chromvar_OPC,genes=rownames(OPC.motif.ctrlvC9ALSFTLD)[1:50],
              show_colnames = F, show_rownames = T,cluster_rows=T,complex=F,
                                      annot.colors=diagnosis_colors, name="z-score",annot.by="diagnosis", 
                                      border_color="lightgray",treeheight_row=20,fontsize=10,show_legend=F,order.by="diagnosis")
ggsave('figures/OPC_chromvar_C9vsCTRL.png', plot=OPC_chromvar_C9vsCTRL, width = 4.5, height = 6, units='in')

OPC_motif<-MotifPlot(
  object = intOPC,
  motifs = head(rownames(OPC.motif.ctrlvC9ALSFTLD),n=10),
  assay = 'peaks')
ggsave('figures/OPC_motif_C9vsCTRL.png',plot=OPC_motif,width=8, height=6, units='in')

astro_chromvar_C9vsCTRL<-dittoHeatmap(avg_chromvar_astro,genes=rownames(astro.motif.ctrlvC9ALSFTLD)[1:50],
                 show_colnames = F, show_rownames = T,cluster_rows=T,complex=F,
                                      annot.colors=diagnosis_colors, name="z-score",annot.by="diagnosis", 
                                      border_color="lightgray",treeheight_row=20,fontsize=10,show_legend=F,order.by="diagnosis")
ggsave('figures/astro_chromvar_C9vsCTRL.png', plot=astro_chromvar_C9vsCTRL, width = 4.5, height = 6, units='in')

astro_motif<-MotifPlot(
  object = intAST,
  motifs = head(rownames(astro.motif.ctrlvC9ALSFTLD),n=10),
  assay = 'peaks')
ggsave('figures/astro_motif_C9vsCTRL.png',plot=astro_motif,width=8, height=6, units='in')

micro_chromvar_C9vsCTRL<-dittoHeatmap(avg_chromvar_micro,genes=rownames(micro.motif.ctrlvC9ALSFTLD)[1:50],
          show_colnames = F, show_rownames = T,cluster_rows=T,complex=F,
                                      annot.colors=diagnosis_colors, name="z-score",annot.by="diagnosis", 
                                      border_color="lightgray",treeheight_row=20,fontsize=10,show_legend=F,order.by="diagnosis")
ggsave('figures/micro_chromvar_C9vsCTRL.png', plot=micro_chromvar_C9vsCTRL, width = 4.5, height = 6, units='in')

micro_motif<-MotifPlot(
  object = intMICRO,
  motifs = head(rownames(micro.motif.ctrlvC9ALSFTLD),n=10),
  assay = 'peaks')
ggsave('figures/micro_motif_C9vsCTRL.png',plot=micro_motif,width=8, height=6, units='in')

rownames(exc.motif.ctrlvC9ALSFTLD)[rownames(exc.motif.ctrlvC9ALSFTLD) == "SMAD2::SMAD3::SMAD4"] <- "SMAD2::4"
exc_chromvar_C9vsCTRL<-dittoHeatmap(avg_chromvar_excitatory,genes=rownames(exc.motif.ctrlvC9ALSFTLD)[1:50],
            show_colnames = F, show_rownames = T,cluster_rows=T,complex=F,
                                      annot.colors=diagnosis_colors, name="z-score",annot.by="diagnosis", 
                                      border_color="lightgray",treeheight_row=20,fontsize=10,show_legend=F,order.by="diagnosis")
ggsave('figures/exc_chromvar_C9vsCTRL.png', plot=exc_chromvar_C9vsCTRL, width = 5, height = 6, units='in')

exc_motif<-MotifPlot(
  object = intEXC,
  motifs = head(rownames(exc.motif.ctrlvC9ALSFTLD),n=10),
  assay = 'peaks')
ggsave('figures/exc_motif_C9vsCTRL.png',plot=exc_motif,width=8, height=6, units='in')

inh_chromvar_C9vsCTRL<-dittoHeatmap(avg_chromvar_inhibitory,genes=rownames(inh.motif.ctrlvC9ALSFTLD)[1:50],
        show_colnames = F, show_rownames = T,cluster_rows=T,complex=F,
                                      annot.colors=diagnosis_colors, name="z-score",annot.by="diagnosis", 
                                      border_color="lightgray",treeheight_row=20,fontsize=10,show_legend=F,order.by="diagnosis")
ggsave('figures/inh_chromvar_C9vsCTRL.png', plot=inh_chromvar_C9vsCTRL, width = 4.5, height = 6, units='in')

inh_motif<-MotifPlot(
  object = intINH,
  motifs = head(rownames(inh.motif.ctrlvC9ALSFTLD),n=10),
  assay = 'peaks')
ggsave('figures/inh_motif_C9vsCTRL.png',plot=inh_motif,width=8, height=6, units='in')

top10<-Extract_Top_Markers(
celltype_motif_markers,
num_genes = 10,
diagnosis_by = "cluster",
rank_by = "avg_log2FC",
gene_column = "gene",
gene_rownames_to_column = FALSE,
data_frame = FALSE,
named_vector = TRUE,
make_unique = FALSE
)

top5<-Extract_Top_Markers(
  celltype_motif_markers,
  num_genes = 5,
  diagnosis_by = "cluster",
  rank_by = "avg_log2FC",
  gene_column = "gene",
  gene_rownames_to_column = FALSE,
  data_frame = FALSE,
  named_vector = TRUE,
  make_unique = FALSE
)


dittoHeatmap(avg_chromvar_ATAC,genes=top5,main = "ATAC DE motifs by celltype",
             annot.by=c("Celltype","diagnoses"),
             show_colnames = F, show_rownames = T,cluster_rows=T,complex=T,
             annot.colors=diagnosis_colors, name="z-score", 
             border_color="lightgray")#,heatmap.colors=colorRampPalette(c("#440154FF","#238A8DFF", "#FDE725FF"))(50))


Idents(microATAC)<-'diagnosis'
C9markers<-FindMarkers(microATAC,ident.1="C9ALSFTLD",ident.2="control")
C9markers
C9markers$gene<-rownames(C9markers)
C9markers
motif_names<-ConvertMotifID(Motifs(intATAC),id=C9markers$gene,assay='peaks')
C9markers$gene<-motif_names
C9markers


oligo_chromvar_C9vsCTRL<-dittoHeatmap(avg_chromvar_oligo,genes=rownames(oligo.motif.ctrlvC9ALSFTLD)[1:50],
                                      show_colnames = F, show_rownames = T,cluster_rows=T,complex=F,
                                      annot.colors=diagnosis_colors, annot.by="diagnosis",name="z-score", 
                                      border_color="lightgray",treeheight_row=20,fontsize=5,fontsize_col=7,show_legend=T,order.by="diagnosis")+NoLegend()

OPC_chromvar_C9vsCTRL<-dittoHeatmap(avg_chromvar_OPC,genes=rownames(OPC.motif.ctrlvC9ALSFTLD)[1:50],
                                      show_colnames = F, show_rownames = T,cluster_rows=T,complex=F,
                                      annot.colors=diagnosis_colors, annot.by="diagnosis",name="z-score", 
                                      border_color="lightgray",treeheight_row=20,fontsize=5,fontsize_col=7,show_legend=T,order.by="diagnosis")+NoLegend()

astro_chromvar_C9vsCTRL<-dittoHeatmap(avg_chromvar_astro,genes=rownames(astro.motif.ctrlvC9ALSFTLD)[1:50],
                                      show_colnames = F, show_rownames = T,cluster_rows=T,complex=F,
                                      annot.colors=diagnosis_colors, annot.by="diagnosis",name="z-score", 
                                      border_color="lightgray",treeheight_row=20,fontsize=5,fontsize_col=7,show_legend=T,order.by="diagnosis")+NoLegend()

micro_chromvar_C9vsCTRL<-dittoHeatmap(avg_chromvar_micro,genes=rownames(micro.motif.ctrlvC9ALSFTLD)[1:50],
                                      show_colnames = F, show_rownames = T,cluster_rows=T,complex=F,
                                      annot.colors=diagnosis_colors, annot.by="diagnosis",name="z-score", 
                                      border_color="lightgray",treeheight_row=20,fontsize=5,fontsize_col=7,show_legend=T,order.by="diagnosis")+NoLegend()

exc_chromvar_C9vsCTRL<-dittoHeatmap(avg_chromvar_exc,genes=rownames(exc.motif.ctrlvC9ALSFTLD)[1:50],
                                      show_colnames = F, show_rownames = T,cluster_rows=T,complex=F,
                                      annot.colors=diagnosis_colors, annot.by="diagnosis",name="z-score", 
                                      border_color="lightgray",treeheight_row=20,fontsize=5,fontsize_col=7,show_legend=T,order.by="diagnosis")+NoLegend()

inh_chromvar_C9vsCTRL<-dittoHeatmap(avg_chromvar_inh,genes=rownames(inh.motif.ctrlvC9ALSFTLD)[1:50],
                                      show_colnames = F, show_rownames = T,cluster_rows=T,complex=F,
                                      annot.colors=diagnosis_colors, annot.by="diagnosis",name="z-score", 
                                      border_color="lightgray",treeheight_row=20,fontsize=5,fontsize_col=7,show_legend=T,order.by="diagnosis")+NoLegend()
