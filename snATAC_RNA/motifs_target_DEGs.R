library(Seurat)
library(Signac)
library(cicero)
library(ggplot2)

intATAC <- readRDS("objects/ATAC_final.RDS")
fragments<-CreateFragmentObject(path="fragments/fragments.tsv.gz",cells=colnames(intATAC))
Fragments(intATAC)<-NULL
Fragments(intATAC)<-fragments
intATAC$diagnosis<-intATAC$group
intATAC$diagnosis<-factor(intATAC$diagnosis,levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))

intOLIGO <- readRDS("objects/ATAC_oligo.RDS")
fragments<-CreateFragmentObject(path="fragments/fragments.tsv.gz",cells=colnames(intOLIGO))
Fragments(intOLIGO)<-NULL
Fragments(intOLIGO)<-fragments
intOLIGO$diagnosis<-intOLIGO$group
intOLIGO$diagnosis<-factor(intOLIGO$diagnosis,levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))

intOPC <- readRDS("objects/ATAC_OPC.RDS")
fragments<-CreateFragmentObject(path="fragments/fragments.tsv.gz",cells=colnames(intOPC))
Fragments(intOPC)<-NULL
Fragments(intOPC)<-fragments
intOPC$diagnosis<-intOPC$group
intOPC$diagnosis<-factor(intOPC$diagnosis,levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))

intAST <- readRDS("objects/ATAC_astrocytes.RDS")
fragments<-CreateFragmentObject(path="fragments/fragments.tsv.gz",cells=colnames(intAST))
Fragments(intAST)<-NULL
Fragments(intAST)<-fragments
intAST$diagnosis<-intAST$group
intAST$diagnosis<-factor(intAST$diagnosis,levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))

intMICRO <- readRDS("objects/ATAC_microglia.RDS")
fragments<-CreateFragmentObject(path="fragments/fragments.tsv.gz",cells=colnames(intMICRO))
Fragments(intMICRO)<-NULL
Fragments(intMICRO)<-fragments
intMICRO$diagnosis<-intMICRO$group
intMICRO$diagnosis<-factor(intMICRO$diagnosis,levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))

intEXC <- readRDS("objects/ATAC_excitatory.RDS")
fragments<-CreateFragmentObject(path="fragments/fragments.tsv.gz",cells=colnames(intEXC))
Fragments(intEXC)<-NULL
Fragments(intEXC)<-fragments
intEXC$diagnosis<-intEXC$group
intEXC$diagnosis<-factor(intEXC$diagnosis,levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))

intINH <- readRDS("objects/ATAC_inhibitory.RDS")
fragments<-CreateFragmentObject(path="fragments/fragments.tsv.gz",cells=colnames(intINH))
Fragments(intINH)<-NULL
Fragments(intINH)<-fragments
intINH$diagnosis<-intINH$group
intINH$diagnosis<-factor(intINH$diagnosis,levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))

#OPC
NFE2L1_motif<-'NFE2L1'
DefaultAssay(intOPC)<-'peaks'
NFE2L1_accessible<-Motifs(intOPC)@data[,NFE2L1_motif]
table(NFE2L1_accessible)
NFE2L1_accessible<-names(NFE2L1_accessible)[NFE2L1_accessible>0]
NFE2L1_accessible_granges<-StringToGRanges(NFE2L1_accessible)
NFE2L1_accessible_granges
NFE2L1_promoters<-promoters(NFE2L1_accessible_granges,upstream=2000,downstream=200)
NFE2L1_promoters
genes_micro<-Annotation(intOPC)
micro_promoters<-promoters(genes_micro,upstream=2000,downstream=200)
micro_promoters
NFE2L1_target_genes<-micro_promoters[micro_promoters@ranges %in% NFE2L1_promoters@ranges]
DE_micro_NFE2L1<-intersect(DE_micro_genes,NFE2L1_target_genes$gene_name)
DE_micro_NFE2L1
opcRNA<-AddModuleScore(opcRNA, features=DE_micro_NFE2L1, name = "NFE2L1_module",ctrl=5)
NFE2L1_module_OPC<-VlnPlot_scCustom(opcRNA,features="NFE2L1_module1",group.by="diagnosis",colors_use=group_colors, pt.size=0)+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("NFE2L1 TF Module")+labs(x="",y="Module z-score")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=6), axis.title.y = element_text(colour="black", size=6))
ggsave('NFE2L1_module_OPC.pdf', plot=NFE2L1_module_OPC, width = 100, height = 80, units='mm')
NFE2L1_module_OPC

#microglia
SPIC_motif<-'SPIC'
DefaultAssay(intMICRO)<-'peaks'
SPIC_accessible<-Motifs(intMICRO)@data[,SPIC_motif]
table(SPIC_accessible)
SPIC_accessible<-names(SPIC_accessible)[SPIC_accessible>0]
SPIC_accessible_granges<-StringToGRanges(SPIC_accessible)
SPIC_accessible_granges
SPIC_promoters<-promoters(SPIC_accessible_granges,upstream=2000,downstream=200)
SPIC_promoters
genes_micro<-Annotation(intMICRO)
micro_promoters<-promoters(genes_micro,upstream=2000,downstream=200)
micro_promoters
SPIC_target_genes<-micro_promoters[micro_promoters@ranges %in% SPIC_promoters@ranges]
DE_micro_SPIC<-intersect(DE_micro_genes,SPIC_target_genes$gene_name)
DE_micro_SPIC
microRNA<-AddModuleScore(microRNA, features=DE_micro_SPIC, name = "SPIC_module",ctrl=5)
SPIC_module_micro<-VlnPlot_scCustom(microRNA,features="SPIC_module1",group.by="diagnosis",colors_use=group_colors, pt.size=0)+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("SPIC TF Module")+labs(x="",y="Module z-score")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=6), axis.title.y = element_text(colour="black", size=6))
ggsave('SPIC_module_micro.pdf', plot=SPIC_module_micro, width = 100, height = 80, units='mm')

CTCF_motif<-'CTCF'
DefaultAssay(intMICRO)<-'peaks'
CTCF_accessible<-Motifs(intMICRO)@data[,CTCF_motif]
table(CTCF_accessible)
CTCF_accessible<-names(CTCF_accessible)[CTCF_accessible>0]
CTCF_accessible_granges<-StringToGRanges(CTCF_accessible)
CTCF_accessible_granges
CTCF_promoters<-promoters(CTCF_accessible_granges,upstream=2000,downstream=200)
CTCF_promoters
genes_micro<-Annotation(intMICRO)
micro_promoters<-promoters(genes,upstream=2000,downstream=200)
micro_promoters
CTCF_target_genes<-micro_promoters[micro_promoters@ranges %in% CTCF_promoters@ranges]
CTCF_target_genes
DE_micro_CTCF<-intersect(DE_micro_genes,CTCF_target_genes$gene_name)
DE_micro_CTCF
microRNA<-AddModuleScore(microRNA, features=DE_micro_CTCF, name = "CTCF_module",ctrl=5)
CTCF_module_micro<-VlnPlot_scCustom(microRNA,features="CTCF_module1",group.by="diagnosis",colors_use=group_colors, pt.size=0)+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("CTCF TF Module")+labs(x="",y="Module z-score")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=6), axis.title.y = element_text(colour="black", size=6))
ggsave('CTCF_module_micro.pdf', plot=CTCF_module_micro, width = 100, height = 80, units='mm')


#astrocytes
FOS_motif<-'FOS'
DefaultAssay(intAST)<-'peaks'
FOS_accessible<-Motifs(intAST)@data[,FOS_motif]
table(FOS_accessible)
FOS_accessible<-names(FOS_accessible)[FOS_accessible>0]
FOS_accessible_granges<-StringToGRanges(FOS_accessible)
FOS_accessible_granges
FOS_promoters<-promoters(FOS_accessible_granges,upstream=2000,downstream=200)
FOS_promoters
genes_AST<-Annotation(intAST)
AST_promoters<-promoters(genes_AST,upstream=2000,downstream=200)
AST_promoters
FOS_target_genes<-AST_promoters[AST_promoters@ranges %in% FOS_promoters@ranges]
FOS_target_genes
DE_ASTPP_FOS<-intersect(DE_ASTPP_genes,FOS_target_genes$gene_name)
DE_ASTPP_FOS
astppRNA<-AddModuleScore(astppRNA, features=DE_ASTPP_FOS, name = "FOS_module",ctrl=5)
FOS_module_ASTPP<-VlnPlot_scCustom(astppRNA,features="FOS_module1",group.by="diagnosis",colors_use=group_colors, pt.size=0)+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("FOS TF Module")+labs(x="",y="Module z-score")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=6), axis.title.y = element_text(colour="black", size=6))
ggsave('FOS_module_ASTPP.pdf', plot=FOS_module_ASTPP, width = 100, height = 80, units='mm')
FOS_module_ASTPP

#EXCITATORY
RELA_motif<-'RELA'
DefaultAssay(intEXC)<-'peaks'
RELA_accessible<-Motifs(intEXC)@data[,RELA_motif]
table(RELA_accessible)
RELA_accessible<-names(RELA_accessible)[RELA_accessible>0]
RELA_accessible_granges<-StringToGRanges(RELA_accessible)
RELA_accessible_granges
RELA_promoters<-promoters(RELA_accessible_granges,upstream=2000,downstream=200)
RELA_promoters
genes_EXC<-Annotation(intEXC)
EXC_promoters<-promoters(genes_EXC,upstream=2000,downstream=200)
EXC_promoters
RELA_target_genes<-EXC_promoters[EXC_promoters@ranges %in% RELA_promoters@ranges]
RELA_target_genes
DE_L23_RELA<-intersect(DE_L23_genes,RELA_target_genes$gene_name)
DE_L23_RELA
L23RNA<-AddModuleScore(L23RNA, features=DE_L23_RELA, name = "RELA_module",ctrl=5)
RELA_module_L23<-VlnPlot_scCustom(L23RNA,features="RELA_module1",group.by="diagnosis",colors_use=group_colors, pt.size=0)+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("RELA TF Module")+labs(x="",y="Module z-score")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=6), axis.title.y = element_text(colour="black", size=6))
ggsave('RELA_module_L23.pdf', plot=RELA_module_L23, width = 100, height = 80, units='mm')
RELA_module_L23

#INHIBITORY
NFYA_motif<-'NFYA'
DefaultAssay(intINH)<-'peaks'
NFYA_accessible<-Motifs(intINH)@data[,NFYA_motif]
table(NFYA_accessible)
NFYA_accessible<-names(NFYA_accessible)[NFYA_accessible>0]
NFYA_accessible_granges<-StringToGRanges(NFYA_accessible)
NFYA_accessible_granges
NFYA_promoters<-promoters(NFYA_accessible_granges,upstream=2000,downstream=200)
NFYA_promoters
genes_INH<-Annotation(intINH)
INH_promoters<-promoters(genes_INH,upstream=2000,downstream=200)
INH_promoters
NFYA_target_genes<-INH_promoters[INH_promoters@ranges %in% NFYA_promoters@ranges]
NFYA_target_genes
DE_INSST_NFYA<-intersect(DE_INSST_genes,NFYA_target_genes$gene_name)
DE_INSST_NFYA
sstRNA<-AddModuleScore(sstRNA, features=DE_INSST_NFYA, name = "NFYA_module",ctrl=5)
NFYA_module_INSST<-VlnPlot_scCustom(sstRNA,features="NFYA_module1",group.by="diagnosis",colors_use=group_colors, pt.size=0)+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("NFYA TF Module")+labs(x="",y="Module z-score")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=6), axis.title.y = element_text(colour="black", size=6))
ggsave('NFYA_module_INSST.pdf', plot=NFYA_module_INSST, width = 100, height = 80, units='mm')
NFYA_module_INSST

DE_INVIP_NFYA<-intersect(DE_INVIP_genes,NFYA_target_genes$gene_name)
DE_INVIP_NFYA
vipRNA<-AddModuleScore(vipRNA, features=DE_INVIP_NFYA, name = "NFYA_module",ctrl=5)
NFYA_module_INVIP<-VlnPlot_scCustom(vipRNA,features="NFYA_module1",group.by="diagnosis",colors_use=group_colors, pt.size=0)+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("NFYA TF Module")+labs(x="",y="Module z-score")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=6), axis.title.y = element_text(colour="black", size=6))
ggsave('NFYA_module_INVIP.pdf', plot=NFYA_module_INVIP, width = 100, height = 80, units='mm')
NFYA_module_INVIP



### old approach


group_colors<-c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2")
celltype_colors<-c("#D55E00","#CC79A7","#E69F00","#0072B2","#009F73","#F0E442")
cellsubtype_colors <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#CC79A7","#666666","#AD7700","#1C91D4","#A04700","#B14380","#4D4D4D","#FFBE2D")

#C9orf72 but includes MOB3B downstream and lots of upstream
CoveragePlot(intATAC,group.by="diagnosis",region="chr9-27325209-27642852",assay='peaks')&scale_fill_manual(values=group_colors)

#C9orf72 with full peak at 3'UTR
CoveragePlot(intATAC,group.by="diagnosis",region="chr9-27545000-27575000",assay='peaks')&scale_fill_manual(values=group_colors)&
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),plot.title = element_text(size=5, face="bold"),
  panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
  legend.title = element_text(colour="black", size=5, face="bold"),
  axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
  axis.title.x = element_text(colour="black", size=5, face="bold"), axis.title.y = element_text(colour="black", size=5, face="bold"))&
  theme(panel.grid.major = element_line(colour="lightgray", size=0.25))

CoveragePlot(intOLIGO,group.by="diagnosis",region="chr9-27545000-27575000",assay='peaks')&scale_fill_manual(values=group_colors)&
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),plot.title = element_text(size=5, face="bold"),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        legend.title = element_text(colour="black", size=5, face="bold"),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=5, face="bold"), axis.title.y = element_text(colour="black", size=5, face="bold"))&
  theme(panel.grid.major = element_line(colour="lightgray", size=0.25))

CoveragePlot(intOPC,group.by="diagnosis",region="chr9-27545000-27575000",assay='peaks')&scale_fill_manual(values=group_colors)&
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),plot.title = element_text(size=5, face="bold"),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        legend.title = element_text(colour="black", size=5, face="bold"),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=5, face="bold"), axis.title.y = element_text(colour="black", size=5, face="bold"))&
  theme(panel.grid.major = element_line(colour="lightgray", size=0.25))

CoveragePlot(intAST,group.by="diagnosis",region="chr9-27545000-27575000",assay='peaks')&scale_fill_manual(values=group_colors)&
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),plot.title = element_text(size=5, face="bold"),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        legend.title = element_text(colour="black", size=5, face="bold"),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=5, face="bold"), axis.title.y = element_text(colour="black", size=5, face="bold"))&
  theme(panel.grid.major = element_line(colour="lightgray", size=0.25))

CoveragePlot(intMICRO,group.by="diagnosis",region="chr9-27545000-27575000",assay='peaks')&scale_fill_manual(values=group_colors)&
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),plot.title = element_text(size=5, face="bold"),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        legend.title = element_text(colour="black", size=5, face="bold"),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=5, face="bold"), axis.title.y = element_text(colour="black", size=5, face="bold"))&
  theme(panel.grid.major = element_line(colour="lightgray", size=0.25))

CoveragePlot(intEXC,group.by="diagnosis",region="chr9-27545000-27575000",assay='peaks')&scale_fill_manual(values=group_colors)&
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),plot.title = element_text(size=5, face="bold"),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        legend.title = element_text(colour="black", size=5, face="bold"),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=5, face="bold"), axis.title.y = element_text(colour="black", size=5, face="bold"))&
  theme(panel.grid.major = element_line(colour="lightgray", size=0.25))

CoveragePlot(intINH,group.by="diagnosis",region="chr9-27545000-27575000",assay='peaks')&scale_fill_manual(values=group_colors)&
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),plot.title = element_text(size=5, face="bold"),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        legend.title = element_text(colour="black", size=5, face="bold"),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=5, face="bold"), axis.title.y = element_text(colour="black", size=5, face="bold"))&
  theme(panel.grid.major = element_line(colour="lightgray", size=0.25))


astro_links <- readRDS("objects/cicero/celltype/intATAC_cicero_links_unadded_astro.RDS")
Links(atac_astro)<-astro_links
atac_astro
DefaultAssay(atac_astro)<-"peaks"
Links(atac_astro)<-astro_links
CoveragePlot(atac_astro,group.by="group",region="chr9-27325209-27642852",assay='peaks')&scale_fill_manual(values=group_colors)
CoveragePlot(atac_astro,group.by="group",region="C9orf72",assay='peaks')&scale_fill_manual(values=group_colors)
intMICRO <- readRDS("objects/ATAC_microglia.RDS")
fragments<-CreateFragmentObject(path="fragments/fragments.tsv.gz",cells=colnames(intMICRO))
Fragments(intATAC)<-NULL
Fragments(intATAC)<-fragments

CoveragePlot(atac_micro,group.by="group",region="C9orf72",assay='peaks')&scale_fill_manual(values=group_colors)
CoveragePlot(atac_micro,group.by="diagnosis",region="C9orf72",assay='peaks')&scale_fill_manual(values=group_colors)
atac_micro$group<-factor(atac_micro$group,levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))
Idents(atac_micro)<-'group'
CoveragePlot(atac_micro,group.by="group",region="C9orf72",assay='peaks')&scale_fill_manual(values=group_colors)
links_micro <- readRDS("objects/cicero/celltype/intATAC_cicero_links_unadded_micro.RDS")
DefaultAssay(atac_micro)<-"peaks"
Links(links_micro)<-links_micro
Links(atac_micro)<-links_micro
CoveragePlot(atac_micro,group.by="group",region="C9orf72",assay='peaks')&scale_fill_manual(values=group_colors)
CoveragePlot(atac_micro,group.by="group",region="TARDBP",assay='peaks')&scale_fill_manual(values=group_colors)
CoveragePlot(atac_micro,group.by="group",region="CD74",assay='peaks')&scale_fill_manual(values=group_colors)
CoveragePlot(atac_micro,group.by="group",region="CSF1R",assay='peaks')&scale_fill_manual(values=group_colors)
micro_links
links_micro
CTCF_motif<-'CTCF'
motif_names<-GetMotifData(atac_micro,slot='motif.names')
motif_names
CTCF_ID<-'MA0139.1'
CTCF_accessible<-Motifs(atac_micro)@data[,CTCF_ID]
atac_micro
DefaultAssay(atac_micro)<-'chromvar'
CTCF_accessible<-Motifs(atac_micro)@data[,CTCF_ID]
DefaultAssay(atac_micro)<-'peaks'
CTCF_accessible<-Motifs(atac_micro)@data[,CTCF_motif]
CTCF_accessible
CTCF_accessible<-names(CTCF_accessible)[CTCF_accessible>0]
CTCF_accessible<-CTCF_accessible[CTCF_accessible %in% VariableFeatures(atac_micro)]
CTCF_accessible
atac_micro
atac_micro<-FindVariableFeatures(atac_micro)
gc()
atac_micro
CTCF_accessible<-Motifs(atac_micro)@data[,CTCF_motif]
CTCF_accessible<-names(CTCF_accessible)[CTCF_accessible>0]
library(Signac)
library(GenomicRanges)
# get gene annotations
genes <- Annotations(atac_micro)
promoters <- promoters(genes, upstream = 2000, downstream = 0)
library(Signac)
library(GenomicRanges)
# get gene annotations
genes <- Annotation(atac_micro)
promoters <- promoters(genes, upstream = 2000, downstream = 0)
promoters
CTCF_target_genes<-CTCF_accessible[CTCF_accessible %in% promoters]
CTCF_accessible
CTCF_target_genes<-promoters$gene_name[match(CTCF_accessible,promoters@ranges)]
head(promoters)
promoters
promoters
CTCF_accessible
promoters
promoters_string<-GRangesToString(promoters)
promoters_string
CTCF_target_genes<-promoters$gene_name[match(CTCF_accessible,promoters_string)]
CTCF_target_genes
table(CTCF_target_genes)
CTCF_target_regions<-match(CTCF_accessible,promoters_string)]
CTCF_target_regions<-match(CTCF_accessible,promoters_string)
CTCF_target_regions
tail(CTCF_target_regions)
CTCF_target_regions<-CTCF_accessible[CTCF_accessible %in% promoters_string, ]
CTCF_accessible
promoters_stringoters
promoters_string
CTCF_accessible
genes<-ClosestFeature(atac_micro, promoters_string)
genes
CTCF_target_regions<-CTCF_accessible[CTCF_accessible %in% promoters_string, ]
CTCF_target_regions<-subset(CTCF_accessible, promoters_string)
CTCF_target_regions<-CTCF_accessible[promoters_string]
CTCF_target_regions
CTCF_target_regions<-CTCF_accessible[promoters_string, ]
promoters_string
CTCF_accessible
identical(CTCF_accessible,promoters_string)
match(CTCF_accessible,promoters_string)
sum(is.na(match(CTCF_accessible,promoters_string)))
promoters_string
summary(CTCF_accessible)
summary(CTCF_motif)
summary(atac_astro
)
intersect(CTCF_accessible,promoters_string)
CTCF_accessible
SPI1_motif<-"SPI1"
SPI1_accessible<-Motifs(atac_micro)@data[,]
SPI1_motif<-"SPI1"
SPI1_accessible<-Motifs(atac_micro)@data[,SPI1_motif]
SPI1_accessible
SPI1_accessible<-names(SPI1_accessible)[SPI1_accessible>0]
SPI1_accessible
promoters
intersect(promoters_string,SPI1_accessible)
promoters
prom_string<-GRangesToString(grange = promoters)
prom_string
intersect(prom_string,SPI1_accessible)
SPI1_promoter<-SPI1_accessible[prom_string]
SPI1_promoter
test<-match(SPI1_accessible,SPI1_promoter)
test
closest <- distanceToNearest(atac_micro, promoters)
closest
promoters <- promoters(genes, upstream = 2000, downstream = 200)
genes
promoters
promoters <- promoters(genes, upstream = 2000, downstream = 200)
library(Signac)
library(GenomicRanges)
# get gene annotations
genes <- Annotation(atac_micro)
promoters <- promoters(genes, upstream = 2000, downstream = 0)
promoters
library(Signac)
library(GenomicRanges)
# get gene annotations
genes <- Annotation(atac_micro)
promoters <- promoters(genes, upstream = 2000, downstream = 200)
test<-match(SPI1_accessible,SPI1_promoter)
test
is.na(sum(test))
is.na(summary(test))
FeaturePlot(atac_astro,features="FOS")
FeaturePlot(atac_astro,features="FOS",assay="chromvar")
DefaultAssay(atac_astro)<-"chromvar"
FeaturePlot(atac_astro,features="FOS")
FeaturePlot(atac_astro,features="FOS",split.by="group")
FeaturePlot_scCustom(atac_astro,features="FOS",split.by="group")
FeaturePlot_scCustom(atac_micro,features="SPIC",split.by="group")
atac_micro
DefaultAssay(atac_micro)<-"chromvar"
FeaturePlot_scCustom(atac_micro,features="SPIC",split.by="group")
FeaturePlot_scCustom(atac_micro,features="CTCF",split.by="group")
FeaturePlot_scCustom(atac_micro,features="SPI1",split.by="group")

