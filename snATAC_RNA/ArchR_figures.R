library(Seurat)
library(Signac)
library(ArchR)
library(ggplot2)
library(patchwork)
library(scCustomize)
library(dittoSeq)
library(ggplotify)
library(ggpubr)
library(viridis)
set.seed(1234)

celltype_colors<-c("#D55E00","#CC79A7","#E69F00","#0072B2","#009F73","#F0E442")
cellsubtype_colors<- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#CC79A7","#666666","#AD7700","#1C91D4","#A04700","#B14380","#4D4D4D","#FFBE2D")

######################
### Signac TSS and Nucleosomal plots
######################

tssplot<-TSSPlot(intATAC,assay="ATAC",group.by="sample")+ggtitle("TSS Enrichment by sample")
ggsave('figures/TSSEnrich_sample.png', plot=tssplot, width = 12, height = 8, dpi=300, units='in')
tssplot<-TSSPlot(intATAC,assay="ATAC",group.by="diagnosis")+ggtitle("TSS Enrichment by diagnosis")
ggsave('figures/TSSEnrich_group.png', plot=tssplot, width = 12, height = 8, dpi=300, units='in')
fragplot<-FragmentHistogram(intATAC,assay="ATAC",group.by="sample")+ggtitle("Nucleosome signal by sample")
ggsave('figures/FragHist_sample.png', plot=fragplot, width = 12, height = 8, dpi=300, units='in')
fragplot<-FragmentHistogram(intATAC,assay="ATAC",group.by="diagnosis")+ggtitle("Nucleosome signal by diagnosis")
ggsave('figures/FragHist_group.png', plot=fragplot, width = 12, height = 8, dpi=300, units='in')

######################
### ArchR TSS and Nucleosomal plots
######################

p1 <- plotGroups(
  ArchRProj = projC9ALS2, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)
p1

p2 <- plotGroups(
  ArchRProj = projC9ALS2, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

p3 <- plotGroups(
  ArchRProj = projC9ALS2, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges"
)

p4 <- plotGroups(
  ArchRProj = projC9ALS2, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

p5 <- plotFragmentSizes(ArchRProj = projC9ALS2)

p6 <- plotTSSEnrichment(ArchRProj = projC9ALS2)

saveRDS(p1, file="snATAC_samples_TSSEnrichment_ridges.RDS")
saveRDS(p2, file="snATAC_samples_TSSEnrichment_violins.RDS")
saveRDS(p3, file="snATAC_samples_log10_nFrags_ridges.RDS")
saveRDS(p4, file="snATAC_samples_log10_nFrags_violins.RDS")
saveRDS(p5, file="snATAC_samples_fragsize_histo.RDS")
saveRDS(p6, file="snATAC_samples_TSSplot.RDS")

ArchR_p1<-p1+scale_color_manual(values=pal_snATAC_ArchR2,guide="none")+
scale_fill_manual(values=pal_snATAC_ArchR2,guide="none")+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
axis.title.x = element_text(colour="black", size=7, face="bold"), axis.title.y = element_text(colour="black", size=7, face="bold"),
axis.text.x = element_text(colour="black", size=7), axis.text.y = element_text(colour="black", size=7))

ArchR_p2<-p2+scale_color_manual(values=pal_snATAC_ArchR,guide="none")+
scale_fill_manual(values=pal_snATAC_ArchR,guide="none")+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
axis.title.x = element_text(colour="black", size=7, face="bold"), axis.title.y = element_text(colour="black", size=7, face="bold"),
axis.text.x = element_text(colour="black", size=7), axis.text.y = element_text(colour="black", size=7))

ArchR_p3<-p3+scale_color_manual(values=pal_snATAC_ArchR2,guide="none")+
scale_fill_manual(values=pal_snATAC_ArchR2,guide="none")+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
axis.title.x = element_text(colour="black", size=7, face="bold"), axis.title.y = element_text(colour="black", size=7, face="bold"),
axis.text.x = element_text(colour="black", size=7), axis.text.y = element_text(colour="black", size=7))

ArchR_p4<-p4+scale_color_manual(values=pal_snATAC_ArchR,guide="none")+
scale_fill_manual(values=pal_snATAC_ArchR,guide="none")+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
axis.title.x = element_text(colour="black", size=7, face="bold"), axis.title.y = element_text(colour="black", size=7, face="bold"),
axis.text.x = element_text(colour="black", size=7), axis.text.y = element_text(colour="black", size=7))

ArchR_p5<-p5+scale_color_manual(values=pal_snATAC_ArchR)+theme(legend.position="none")+
xlab("snATAC-Seq Fragment Size (bp)")+
ylab("Proportion of Fragments")+
guides(color=guide_legend(ncol=1,title = "Sample",face="bold",size=7))+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
axis.text.x = element_text(colour="black", size=7), axis.text.y = element_text(colour="black", size=7),
axis.title.x = element_text(colour="black", size=7, face="bold"), axis.title.y = element_text(colour="black", size=7, face="bold"))

ArchR_p6<-p6+scale_color_manual(values=pal_snATAC_ArchR)+theme(legend.position="right", legend.text=element_text(size=5), legend.key.size = unit(0, 'lines'))+
xlab("snATAC-Seq Fragment Size (bp)")+
ylab("Proportion of Fragments")+
guides(color=guide_legend(ncol=1,title = "sample"))+
theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
legend.title = element_text(colour="black", size=7, face="bold"),
axis.text.x = element_text(colour="black", size=7), axis.text.y = element_text(colour="black", size=7),
axis.title.x = element_text(colour="black", size=7, face="bold"), axis.title.y = element_text(colour="black", size=7, face="bold"))

ArchR_plots <-wrap_plots(ArchR_p2,ArchR_p4,ArchR_p5,ArchR_p6)+plot_layout(ncol=2,nrow=2)&theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))
ggsave(plot=ArchR_plots, file = "ArchR_plots_FigS1.pdf", width = 180, height = 120, units="mm")

#ArchR_plots <-wrap_plots(ArchR_p2,ArchR_p4,ArchR_p5,ArchR_p6)+plot_layout(ncol=2,nrow=2)+plot_annotation(tag_levels = 'a') &theme(plot.tag = element_text(size=10,face = 'bold'), plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"))
#ggsave(plot=ArchR_plots, file = "ArchR_plots_FigS1.pdf", width = 180, height = 120, units="mm")
