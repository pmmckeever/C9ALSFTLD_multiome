library(Seurat)
library(readr)
library(dplyr)
library(dittoSeq)
set.seed(1234)

intRNA <- readRDS(file="objects/intRNA_V3only.RDS")
Idents(intRNA)<-"cellsubtype"

intL23<-subset(intRNA,idents="L23")
intL4<-subset(intRNA,idents="L4")
intL56<-subset(intRNA,idents="L56")
intL56CC<-subset(intRNA,idents="L56CC")

intINPV<-subset(intRNA,idents="IN-PV")
intINSST<-subset(intRNA,idents="IN-SST")
intINSV2C<-subset(intRNA,idents="IN-SV2C")
intINVIP<-subset(intRNA,idents="IN-VIP")

Xiao2011_targets<-read_csv("ref_data/TDP43_targets.csv")
MA2022_genes <- read_csv("ref_data/MA_2022.csv")
MA2022_genes

intRNA.L23.avg<-AverageExpression(intL23,assays="RNA",return.seurat=T)
intRNA.L4.avg<-AverageExpression(intL4,assays="RNA",return.seurat=T)
intRNA.L56.avg<-AverageExpression(intRNA.ex3,assays="RNA",return.seurat=T)
intRNA.L56CC.avg<-AverageExpression(intRNA.ex4,assays="RNA",return.seurat=T)

intRNA.INPV.avg<-AverageExpression(intPV,assays="RNA",return.seurat=T)
intRNA.INSST.avg<-AverageExpression(intSST,assays="RNA",return.seurat=T)
intRNA.INSV2C.avg<-AverageExpression(intSV2C,assays="RNA",return.seurat=T)
intRNA.INVIP.avg<-AverageExpression(intVIP,assays="RNA",return.seurat=T)

intRNA.ASTPP.avg$diagnosis<-Idents(intRNA.ASTPP.avg)
intRNA.ASTFB.avg$diagnosis<-Idents(intRNA.ASTFB.avg)

intRNA.L23.avg$diagnosis<-Idents(intRNA.L23.avg)
intRNA.L4.avg$diagnosis<-Idents(intRNA.L4.avg)
#intRNA.L56.avg$diagnosis<-Idents(intRNA.L56.avg)
#intRNA.L56CC.avg$diagnosis<-Idents(intRNA.L56CC.avg)

intRNA.INPV.avg$diagnosis<-Idents(intRNA.INPV.avg)
intRNA.INSST.avg$diagnosis<-Idents(intRNA.INSST.avg)
intRNA.INSV2C.avg$diagnosis<-Idents(intRNA.INSV2C.avg)
intRNA.INVIP.avg$diagnosis<-Idents(intRNA.INVIP.avg)

Idents(intL23)<-'diagnosis'
Idents(intL4)<-'diagnosis'
Idents(intPV)<-'diagnosis'
Idents(intSST)<-'diagnosis'
Idents(intSV2C)<-'diagnosis'
Idents(intVIP)<-'diagnosis'

DE_TDP_L23<-cellsubtype.C9ALSFTLD.cutoff$gene[cellsubtype.C9ALSFTLD.cutoff$cellsubtype == "L23"]
DE_TDP_L4<-cellsubtype.C9ALSFTLD.cutoff$gene[cellsubtype.C9ALSFTLD.cutoff$cellsubtype == "L4"]
DE_TDP_L56<-cellsubtype.C9ALSFTLD.cutoff$gene[cellsubtype.C9ALSFTLD.cutoff$cellsubtype == "L56"]
DE_TDP_L56CC<-cellsubtype.C9ALSFTLD.cutoff$gene[cellsubtype.C9ALSFTLD.cutoff$cellsubtype == "L56-CC"]

DE_TDP_INPV<-cellsubtype.C9ALSFTLD.cutoff$gene[cellsubtype.C9ALSFTLD.cutoff$cellsubtype == "IN-PV"]
DE_TDP_INSST<-cellsubtype.C9ALSFTLD.cutoff$gene[cellsubtype.C9ALSFTLD.cutoff$cellsubtype == "IN-SST"]
DE_TDP_INSV2C<-cellsubtype.C9ALSFTLD.cutoff$gene[cellsubtype.C9ALSFTLD.cutoff$cellsubtype == "IN-SV2C"]
DE_TDP_INVIP<-cellsubtype.C9ALSFTLD.cutoff$gene[cellsubtype.C9ALSFTLD.cutoff$cellsubtype == "IN-VIP"]

overlap_L23<-intersect(DE_TDP_L23,pull(Xiao2011_targets))
overlap_L4<-intersect(DE_TDP_L4,pull(Xiao2011_targets))
overlap_L56<-intersect(DE_TDP_L56,pull(Xiao2011_targets))
overlap_L56CC<-intersect(DE_TDP_L56CC,pull(Xiao2011_targets))

overlap_INPV<-intersect(DE_TDP_INPV,pull(Xiao2011_targets))
overlap_INSST<-intersect(DE_TDP_INSST,pull(Xiao2011_targets))
overlap_INSV2C<-intersect(DE_TDP_INSV2C,pull(Xiao2011_targets))
overlap_INVIP<-intersect(DE_TDP_INVIP,pull(Xiao2011_targets))
overlap

overlap2_L23<-intersect(DE_TDP_L23,pull(MA2022_genes))
overlap2_L4<-intersect(DE_TDP_L4,pull(MA2022_genes))
overlap2_L56<-intersect(DE_TDP_L56,pull(MA2022_genes))
overlap2_L56CC<-intersect(DE_TDP_L56CC,pull(MA2022_genes))
overlap2_INPV<-intersect(DE_TDP_INPV,pull(MA2022_genes))
overlap2_INSST<-intersect(DE_TDP_INSST,pull(MA2022_genes))
overlap2_INSV2C<-intersect(DE_TDP_INSV2C,pull(MA2022_genes))
overlap2_INVIP<-intersect(DE_TDP_INVIP,pull(MA2022_genes))

intL23<-AddModuleScore(intL23,features=overlap2_L23,ctrl=5,name="TDPsplice")
intL4<-AddModuleScore(intL4,features=overlap2_L4,ctrl=5,name="TDPsplice")
intPV<-AddModuleScore(intPV,features=overlap2_INPV,ctrl=5,name="TDPsplice")
intSST<-AddModuleScore(intSST,features=overlap2_INSST,ctrl=5,name="TDPsplice")
intSV2C<-AddModuleScore(intSV2C,features=overlap2_INSV2C,ctrl=5,name="TDPsplice")
intVIP<-AddModuleScore(intVIP,features=overlap2_INVIP,ctrl=5,name="TDPsplice")


excite_TDP<-VlnPlot_scCustom(intEXC,features="TDPgenes1",group.by = "group",pt.size=0,colors_use = group_colors)+
geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("DE TDP43 targets")+labs(x="",y="Module z-score")+
theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
axis.title.x = element_text(colour="black", size=5), axis.title.y = element_text(colour="black", size=5))
excite_TDP
ggsave('excite_TDPtargets.pdf', plot=excite_TDP, width = 100, height = 120, units='mm')

#low experssion commented out
intL23<-AddModuleScore(intL23,features=overlap,ctrl=5,name="TDPgenes")
intL4<-AddModuleScore(intL4,features=overlap,ctrl=5,name="TDPgenes")
#intL56<-AddModuleScore(intL56,features=overlap,ctrl=5,name="TDPgenes")
#intL56CC<-AddModuleScore(intL56CC,features=overlap,ctrl=5,name="TDPgenes")

intINPV<-AddModuleScore(intINPV,features=overlap_INPV,ctrl=5,name="TDPgenes")
intINSST<-AddModuleScore(intINSST,features=overlap_INSST,ctrl=5,name="TDPgenes")
#intINSV2C<-AddModuleScore(intINSV2C,features=overlap_INSV2C,ctrl=5,name="TDPgenes")
intINVIP<-AddModuleScore(intINVIP,features=overlap_INVIP,ctrl=5,name="TDPgenes")

excite_L23<-VlnPlot_scCustom(intL23,features="TDPgenes1",group.by = "group",pt.size=0,colors_use = group_colors)+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("L2/3 TDP43 targets")+labs(x="",y="Module z-score")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=5), axis.title.y = element_text(colour="black", size=5))
ggsave('L23_TDPtargets.pdf', plot=excite_L23, width = 100, height = 120, units='mm')

excite_L4<-VlnPlot_scCustom(intL4,features="TDPgenes1",group.by = "group",pt.size=0,colors_use = group_colors)+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("L4 TDP43 targets")+labs(x="",y="Module z-score")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=5), axis.title.y = element_text(colour="black", size=5))
ggsave('L4_TDPtargets.pdf', plot=excite_L4, width = 100, height = 120, units='mm')

inhibit_INPV<-VlnPlot_scCustom(intINPV,features="TDPgenes1",group.by = "group",pt.size=0,colors_use = group_colors)+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("IN-PV TDP43 targets")+labs(x="",y="Module z-score")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
  axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
  axis.title.x = element_text(colour="black", size=5), axis.title.y = element_text(colour="black", size=5))
ggsave('INPV_TDPtargets.pdf', plot=inhibit_INPV, width = 100, height = 120, units='mm')

inhibit_INSST<-VlnPlot_scCustom(intINSST,features="TDPgenes1",group.by = "group",pt.size=0,colors_use = group_colors)+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("IN-SST TDP43 targets")+labs(x="",y="Module z-score")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
  axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
  axis.title.x = element_text(colour="black", size=5), axis.title.y = element_text(colour="black", size=5))
ggsave('INSST_TDPtargets.pdf', plot=inhibit_INSST, width = 100, height = 120, units='mm')

inhibit_INVIP<-VlnPlot_scCustom(intINVIP,features="TDPgenes1",group.by = "group",pt.size=0,colors_use = group_colors)+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("IN-VIP TDP43 targets")+labs(x="",y="Module z-score")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
  axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
  axis.title.x = element_text(colour="black", size=5), axis.title.y = element_text(colour="black", size=5))
ggsave('INVIP_TDPtargets.pdf', plot=inhibit_INVIP, width = 100, height = 120, units='mm')

FOS_astro<-VlnPlot_scCustom(atac_astro,features="FOS",group.by = "group",pt.size=0,colors_use = group_colors)+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("FOS")+labs(x="",y="chromVAR deviations")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=5), axis.title.y = element_text(colour="black", size=5))
ggsave('FOS_astro.pdf', plot=FOS_astro, width = 100, height = 120, units='mm')


NFE2L1_OPC<-VlnPlot_scCustom(intOPC,features="NFE2L1",group.by = "group",pt.size=0,colors_use = group_colors, assay="chromvar")+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("NFE2L1")+labs(x="",y="chromVAR deviations")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=5), axis.title.y = element_text(colour="black", size=5))
ggsave('NFE2L1_OPC.pdf', plot=NFE2L1_OPC, width = 100, height = 120, units='mm')

up_in_DAA<-list(c("CTSB","CSMD1","C4B","VIM","SPARCL1","AQP4","GFAP"))
down_in_DAA<-list(c("SLC7A10","TRPM3","LUZP2"))

intASTPP<-AddModuleScore(intASTPP, features=up_in_DAA, name = "up_in_DAA",ctrl=5)
intASTPP<-AddModuleScore(intASTPP, features=down_in_DAA, name = "down_in_DAA",ctrl=5)

intASTFB<-AddModuleScore(intASTFB, features=up_in_DAA, name = "up_in_DAA",ctrl=5)
intASTFB<-AddModuleScore(intASTFB, features=down_in_DAA, name = "down_in_DAA",ctrl=5)

DAA_astroPP<-VlnPlot_scCustom(intASTPP,features="up_in_DAA1",group.by = "group",pt.size=0,colors_use = group_colors)+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("DAA markers")+labs(x="",y="Module z-score")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=5), axis.title.y = element_text(colour="black", size=5))
ggsave('DAA_astroPP.pdf', plot=DAA_astroPP, width = 100, height = 120, units='mm')

DAA_down_astroPP<-VlnPlot_scCustom(intASTPP,features="down_in_DAA1",group.by = "group",pt.size=0,colors_use = group_colors)+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("Down in DAA")+labs(x="",y="Module z-score")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=5), axis.title.y = element_text(colour="black", size=5))
ggsave('DAA_down_astroPP.pdf', plot=DAA_down_astroPP, width = 100, height = 120, units='mm')

DAA_astroFB<-VlnPlot_scCustom(intASTFB,features="up_in_DAA1",group.by = "group",pt.size=0,colors_use = group_colors)+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("DAA markers")+labs(x="",y="Module z-score")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=5), axis.title.y = element_text(colour="black", size=5))
ggsave('DAA_astroFB.pdf', plot=DAA_astroFB, width = 100, height = 120, units='mm')

DAA_down_astroFB<-VlnPlot_scCustom(intASTFB,features="down_in_DAA1",group.by = "group",pt.size=0,colors_use = group_colors)+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("Down in DAA")+labs(x="",y="Module z-score")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=5), axis.title.y = element_text(colour="black", size=5))
ggsave('DAA_down_astroFB.pdf', plot=DAA_down_astroFB, width = 100, height = 120, units='mm')


####
GFAP_astroPP<-VlnPlot_scCustom(intASTPP,features="GFAP",group.by = "group",pt.size=0,colors_use = group_colors)+
        geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("GFAP")+labs(x="",y="Expression level")+
        theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=9), axis.text.y = element_text(colour="black", size=9),
        axis.title.x = element_text(colour="black", size=10, face="bold"), axis.title.y = element_text(colour="black", size=10, face="bold"))
ggsave('GFAP_astroPP.pdf', plot=GFAP_astroPP, width = 100, height = 120, units='mm')

SLC1A2_astroPP<-VlnPlot_scCustom(intASTPP,features="SLC1A2",group.by = "group",pt.size=0,colors_use = group_colors)+
        geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("SLC1A2")+labs(x="",y="Expression level")+
        theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=9), axis.text.y = element_text(colour="black", size=9),
        axis.title.x = element_text(colour="black", size=10, face="bold"), axis.title.y = element_text(colour="black", size=10, face="bold"))
ggsave('SLC1A2_astroPP.pdf', plot=SLC1A2_astroPP, width = 100, height = 120, units='mm')


AQP4_astroPP<-VlnPlot_scCustom(intASTPP,features="AQP4",group.by = "group",pt.size=0,colors_use = group_colors)+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("AQP4")+labs(x="",y="Expression level")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=9), axis.text.y = element_text(colour="black", size=9),
        axis.title.x = element_text(colour="black", size=10, face="bold"), axis.title.y = element_text(colour="black", size=10, face="bold"))
ggsave('AQP4_astroPP.pdf', plot=AQP4_astroPP, width = 100, height = 120, units='mm')


### TDP splice targets from Ma et al., 2022

L23_TDPsplice<-VlnPlot_scCustom(intL23,features="TDPsplice1",group.by = "group",pt.size=0,colors_use = group_colors)+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("L2/3 TDP43 splice")+labs(x="",y="Module z-score")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=5), axis.title.y = element_text(colour="black", size=5))
L23_TDPsplice
ggsave('L23_TDPsplice.pdf', plot=L23_TDPsplice, width = 100, height = 120, units='mm')

L4_TDPsplice<-VlnPlot_scCustom(intL4,features="TDPsplice1",group.by = "group",pt.size=0,colors_use = group_colors)+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("L4 TDP43 splice")+labs(x="",y="Module z-score")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=5), axis.title.y = element_text(colour="black", size=5))
L4_TDPsplice
ggsave('L4_TDPsplice.pdf', plot=L4_TDPsplice, width = 100, height = 120, units='mm')

INPV_TDPsplice<-VlnPlot_scCustom(intPV,features="TDPsplice1",group.by = "group",pt.size=0,colors_use = group_colors)+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("IN-PV TDP43 splice")+labs(x="",y="Module z-score")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=5), axis.title.y = element_text(colour="black", size=5))
INPV_TDPsplice
ggsave('INPV_TDPsplice.pdf', plot=INPV_TDPsplice, width = 100, height = 120, units='mm')

INSST_TDPsplice<-VlnPlot_scCustom(intSST,features="TDPsplice1",group.by = "group",pt.size=0,colors_use = group_colors)+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("IN-SST TDP43 splice")+labs(x="",y="Module z-score")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=5), axis.title.y = element_text(colour="black", size=5))
INSST_TDPsplice
ggsave('INSST_TDPsplice.pdf', plot=INSST_TDPsplice, width = 100, height = 120, units='mm')

INVIP_TDPsplice<-VlnPlot_scCustom(intVIP,features="TDPsplice1",group.by = "group",pt.size=0,colors_use = group_colors)+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("IN-VIP TDP43 splice")+labs(x="",y="Module z-score")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=5), axis.title.y = element_text(colour="black", size=5))
INVIP_TDPsplice
ggsave('INVIP_TDPsplice.pdf', plot=INVIP_TDPsplice, width = 100, height = 120, units='mm')