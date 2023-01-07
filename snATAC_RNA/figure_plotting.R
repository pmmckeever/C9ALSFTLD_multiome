library(Seurat)
library(Signac)
library(ggplot2)
set.seed(1234)


##################
# AGE BARPLOTS
##################
#snRNA age and condition data
snRNA_age<-c(50,59,72,72,53,48,59,57,65,88,43,72,37,50,60,60,58,47,59,72,71,66,88,89)
snRNA_sex<-c("female","male","female","male","female","female","male","female","male","female","male","female","male","male","male","male","male","female","female","male","female","female","female","male")
snRNA_samples<-c("CTRL1","CTRL2","CTRL3","CTRL4","CTRL5","CTRL6",
"sALSnoFTLD1","sALSnoFTLD2","sALSnoFTLD3","sALSnoFTLD4","sALSnoFTLD5","sALSnoFTLD6","sALSnoFTLD7","sALSnoFTLD8",
"C9ALSFTLD1","C9ALSFTLD2","C9ALSFTLD3","C9ALSFTLD4","C9ALSFTLD5","C9ALSFTLD6",
"C9ALSnoFTLD1","C9ALSnoFTLD2","C9ALSnoFTLD3","C9noALSnoFTLD")
snRNA_diagnosis<-c("control","control","control","control","control","control","sALSnoFTLD","sALSnoFTLD","sALSnoFTLD","sALSnoFTLD","sALSnoFTLD",
"sALSnoFTLD","sALSnoFTLD","sALSnoFTLD","C9ALSFTLD","C9ALSFTLD","C9ALSFTLD","C9ALSFTLD","C9ALSFTLD","C9ALSFTLD","C9ALSnoFTLD",
"C9ALSnoFTLD","C9ALSnoFTLD","C9noALSnoFTLD")
group_colors<-c("control" = "#E69F00","C9noALSnoFTLD"="#56B4E9","C9ALSFTLD"="#009E73","C9ALSnoFTLD"="#F0E442","sALSnoFTLD"="#0072B2")
snRNA_age_data<-data.frame(snRNA_samples,snRNA_diagnosis,snRNA_age,snRNA_sex)

#snATAC age and condition data
snATAC_age<-c(50,59,72,72,59,57,65,88,43,72,37,50,60,60,58,59,72,57,71,66,88,89)
snATAC_sex<-c("female","male","female","male","male","female","male","female","male","female","male","male","male","male","male","female","male","male","female","female","female","male")
snATAC_samples<-c("CTRL1","CTRL2","CTRL3","CTRL4",
"sALSnoFTLD1","sALSnoFTLD2","sALSnoFTLD3","sALSnoFTLD4","sALSnoFTLD5","sALSnoFTLD6","sALSnoFTLD7","sALSnoFTLD8",
"C9ALSFTLD1","C9ALSFTLD2","C9ALSFTLD3","C9ALSFTLD5","C9ALSFTLD6","C9ALSFTLD7",
"C9ALSnoFTLD1","C9ALSnoFTLD2","C9ALSnoFTLD3","C9noALSnoFTLD")
snATAC_diagnosis<-c("control","control","control","control","sALSnoFTLD","sALSnoFTLD","sALSnoFTLD","sALSnoFTLD","sALSnoFTLD",
"sALSnoFTLD","sALSnoFTLD","sALSnoFTLD","C9ALSFTLD","C9ALSFTLD","C9ALSFTLD","C9ALSFTLD","C9ALSFTLD","C9ALSFTLD","C9ALSnoFTLD",
"C9ALSnoFTLD","C9ALSnoFTLD","C9noALSnoFTLD")
group_colors<-c("control" = "#E69F00","C9noALSnoFTLD"="#56B4E9","C9ALSFTLD"="#009E73","C9ALSnoFTLD"="#F0E442","sALSnoFTLD"="#0072B2")
snATAC_age_data<-data.frame(snATAC_samples,snATAC_diagnosis,snATAC_age,snATAC_sex)

snRNA_age_data$tech<-"snRNA"
snATAC_age_data$tech<-"snATAC"

snRNA_age_data$samples<-snRNA_age_data$snRNA_samples
snRNA_age_data$diagnosis<-snRNA_age_data$snRNA_diagnosis
snRNA_age_data$age<-snRNA_age_data$snRNA_age
snRNA_age_data$sex<-snRNA_age_data$snRNA_sex
snATAC_age_data$samples<-snATAC_age_data$snATAC_samples
snATAC_age_data$diagnosis<-snATAC_age_data$snATAC_diagnosis
snATAC_age_data$age<-snATAC_age_data$snATAC_age
snATAC_age_data$sex<-snATAC_age_data$snATAC_sex

snRNA_age_data$snRNA_samples<-NULL
snRNA_age_data$snRNA_diagnosis<-NULL
snRNA_age_data$snRNA_age<-NULL
snRNA_age_data$snRNA_sex<-NULL
snATAC_age_data$snATAC_samples<-NULL
snATAC_age_data$snATAC_diagnosis<-NULL
snATAC_age_data$snATAC_age<-NULL
snATAC_age_data$snATAC_sex<-NULL

age_data<-rbind(snRNA_age_data,snATAC_age_data)

age_data$tech<-factor(age_data$tech,levels=c("snRNA","snATAC"))

#plot both tech ages together

age_plot<-ggplot(age_data, aes(x=tech,y=age))+
geom_boxplot() +geom_jitter(aes(colour=diagnosis,shape=sex), size=4)+
scale_color_manual(values=group_colors,name="Diagnosis")+scale_shape_manual(values=c(19,17),name="Sex")+
labs(x="Samples",y="Age")+theme(legend.position="right",plot.title = element_text(size=7, face="bold"),
          legend.text=element_text(size=7),
          axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
          axis.title.x = element_text(colour="black", size=7, face="bold"), 
          axis.title.y = element_text(colour="black", size=7, face="bold"),
          axis.text.x = element_text(colour="black", size=7), 
          axis.text.y = element_text(colour="black", size=7))+theme(text=element_text(family="Arial"))+
ggtitle("snRNA by Age")+
guides(colour = guide_legend(override.aes = list(shape = 15, size = 4)))

# individual barplots

snRNA_age<-ggplot(snRNA_age_data, aes(x="",y=snRNA_age))+
geom_boxplot() +geom_jitter(aes(colour=snRNA_diagnosis,shape=snRNA_sex), size=5)+
scale_color_manual(values=group_colors,name="Diagnosis")+scale_shape_manual(values=c(19,17),name="Sex")+
labs(x="Samples",y="Age")+theme(legend.position="right",plot.title = element_text(size=7, face="bold"),
          legend.text=element_text(size=7),
          axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
          axis.title.x = element_text(colour="black", size=7, face="bold"), 
          axis.title.y = element_text(colour="black", size=7, face="bold"),
          axis.text.x = element_text(colour="black", size=7), 
          axis.text.y = element_text(colour="black", size=7))+
ggtitle("snRNA by Age")+NoLegend()
#guides(colour = guide_legend(override.aes = list(shape = 15)))
#ggsave('figures/snRNA_age.png', width=4.5, height=8, units='in')
saveRDS(snRNA_age,file="snRNA-Seq_sample_age.RDS")

snATAC_age<-ggplot(snATAC_age_data, aes(x="",y=snATAC_age)) +
geom_boxplot() +geom_jitter(aes(colour=snATAC_diagnosis, shape=snATAC_sex), size=5)+
scale_color_manual(values=group_colors,name="Diagnosis")+scale_shape_manual(values=c(19,17),name="Sex")+
labs(x="Samples",y="Age")+theme(legend.position="right",plot.title = element_text(size=7, face="bold"),
          legend.text=element_text(size=7),
          axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
          axis.title.x = element_text(colour="black", size=7, face="bold"), 
          axis.title.y = element_text(colour="black", size=7, face="bold"),
          axis.text.x = element_text(colour="black", size=7), 
          axis.text.y = element_text(colour="black", size=7))+ggtitle("snRNA by Age")+
guides(colour = guide_legend(override.aes = list(shape = 15, size = 4)))
#ggsave('figures/snATAC_age.png', width=4.5, height=8, units='in', dpi=300)
saveRDS(snATAC_age,file="snATAC-Seq_sample_age.RDS")

intRNA <- readRDS(file="objects/intRNA_final.RDS")
intATAC <- readRDS(file="objects/intATAC_final.RDS")

blank_theme <- theme(
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  plot.background=element_blank()
)

p1<-DimPlot_scCustom(intRNA,group.by="group",pt.size=0.001, colors_use=group_colors)+blank_theme+theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))+theme(title=element_text(size=7))+ggtitle("snRNA Diagnosis")
p2<-DimPlot_scCustom(intATAC,group.by="group",pt.size=0.001, colors_use=group_colors)+blank_theme+theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))+theme(title=element_text(size=7))+ggtitle("snATAC Diagnosis")
pdf(file="group_dimplots.pdf",width=3.6,height =1.8)
p1 + p2
dev.off()

p3<-DimPlot_scCustom(intRNA,group.by="seurat_clusters",pt.size=0.001)+blank_theme+theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))+
theme(title=element_text(size=7), legend.position="right",legend.text=element_text(size=5), legend.key.size = unit(0.5, 'lines'))+
guides(color = guide_legend(override.aes = list(size=2), ncol=2))

p4<-DimPlot_scCustom(intATAC,group.by="seurat_clusters",pt.size=0.001)+blank_theme+theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))+
theme(title=element_text(size=7), legend.position="right",legend.text=element_text(size=5), legend.key.size = unit(0.5, 'lines'))+
guides(color = guide_legend(override.aes = list(size=2), ncol=2))

seurat_clusters_UMAPs<-p3+p4
ggsave("seurat_clusters_UMAPs.tiff", seurat_clusters_UMAPs, width=160, height=80, units='mm', dpi=300)
ggsave("seurat_clusters_UMAPs.pdf", seurat_clusters_UMAPs, width=160, height=80,units='mm')
dev.off()

p5<-DimPlot_scCustom(intRNA,group.by="celltype_clusters",pt.size=0.001,label=T)+blank_theme+theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))+theme(title=element_text(size=7))+NoLegend()
p6<-DimPlot_scCustom(intATAC,group.by="celltype_clusters",pt.size=0.001,label=T)+blank_theme+theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))+theme(title=element_text(size=7))+NoLegend()

pal_snRNA<-c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7",
             "#666666","#AD7700","#1C91D4","#007756","#D5C711","#005685","#B14380","#4D4D4D",
             "#FFBE2D","#80C7EF","#00F6B3","#F4EB71","#06A5FF","#FF8320","#D99BBD","#8C8C8C","#FFCB57")
#22 samples
pal_snATAC<-c("#E69F00","#56B4E9","#009E73","#F0E442","#CC79A7","#666666",
              "#AD7700","#1C91D4","#D5C711","#005685","#A04700","#B14380","#4D4D4D","#FFBE2D",
              "#80C7EF","#00F6B3","#F4EB71","#06A5FF","#FF8320","#D99BBD","#8C8C8C","#FFCB57")

p7<-DimPlot_scCustom(intRNA,group.by="sample",pt.size=0.001,colors_use=pal_snRNA)+blank_theme+theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))+theme(title=element_text(size=7))+NoLegend()
p8<-DimPlot_scCustom(intRNA_batch,group.by="sample",pt.size=0.001, colors_use = pal_snRNA)+blank_theme+theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))+theme(title=element_text(size=7))+NoLegend()
intRNA_batch_UMAPs<-p7+p8
ggsave("intRNA_batch_UMAPs.tiff", intRNA_batch_UMAPs, width=160, height=80, units='mm', dpi=300)
#ggsave("intRNA_batch_UMAPs.pdf", intRNA_batch_UMAPs, width=160, height=80,units='mm')

p9<-DimPlot_scCustom(intATAC,group.by="sample",pt.size=0.001, colors_use = pal_snATAC)+blank_theme+theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))+theme(title=element_text(size=7))+NoLegend()
p10<-DimPlot_scCustom(intATAC_batch,group.by="sample",pt.size=0.001, colors_use = pal_snATAC)+blank_theme+theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))+theme(title=element_text(size=7))+NoLegend()

intATAC_batch_UMAPs<-p9+p10
ggsave("intATAC_batch_UMAPs.tiff", intATAC_batch_UMAPs, width=160, height=80, units='mm', dpi=300)
#ggsave("intATAC_batch_UMAPs.pdf", intATAC_batch_UMAPs, width=160, height=80,units='mm')

p7<-DimPlot_scCustom(intRNA,group.by="sample",pt.size=0.001,colors_use=pal_snRNA)+blank_theme+theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))+theme(title=element_text(size=7))+NoLegend()
p8<-DimPlot_scCustom(intRNA_batch,group.by="sample",pt.size=0.001,colors_use=pal_snRNA)+blank_theme+theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))+
theme(title=element_text(size=7), legend.position="right",legend.text=element_text(size=5), legend.key.size = unit(0.5, 'lines'))+
guides(color = guide_legend(override.aes = list(size=2), ncol=1))

intRNA_batch_UMAPs<-p7+p8
#ggsave("intRNA_batch_UMAPs.tiff", intRNA_batch_UMAPs, width=160, height=80, units='mm', dpi=300)
ggsave("intRNA_batch_UMAPs.pdf", intRNA_batch_UMAPs, width=160, height=80,units='mm')

p9<-DimPlot_scCustom(intATAC,group.by="sample",pt.size=0.001, colors_use = pal_snATAC)+blank_theme+theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))+theme(title=element_text(size=7))+NoLegend()
p10<-DimPlot_scCustom(intATAC_batch,group.by="sample",pt.size=0.001, colors_use = pal_snATAC)+blank_theme+theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))+
theme(title=element_text(size=7), legend.position="right",legend.text=element_text(size=5), legend.key.size = unit(0.5, 'lines'))+
guides(color = guide_legend(override.aes = list(size=2), ncol=1))

intATAC_batch_UMAPs<-p9+p10
#ggsave("intATAC_batch_UMAPs.tiff", intATAC_batch_UMAPs, width=160, height=80, units='mm', dpi=300)
ggsave("intATAC_batch_UMAPs.pdf", intATAC_batch_UMAPs, width=160, height=80,units='mm')


#Sample celltype frequency barplots

sample1_BarPlot_RNA<-dittoBarPlot(intRNA, "cellsubtype",group.by="sample", 
                                       x.reorder=c(11,12,13,14,15,16,10,1,2,3,4,5,6,7,8,9,17,18,19,20,21,22,23,24),
                                       var.labels.reorder=c(12,13,2,1,11,7,8,9,10,3,4,5,6),
                                       color.panel=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7",
                                                     "#666666","#AD7700","#1C91D4","#D5C711","#A04700","#B14380","#4D4D4D","#FFBE2D"),
                                       scale="count")+ggtitle("snRNA sample cell numbers")+labs(x="Sample",y="Number of cells")+
  mytheme+NoLegend()
ggsave('sample1_BarPlot_RNA.pdf', plot=sample1_BarPlot_RNA, width = 6, height = , units='mm')


sample2_BarPlot_RNA<-dittoBarPlot(intRNA, "cellsubtype",group.by="sample", 
                                       x.reorder=c(11,12,13,14,15,16,10,1,2,3,4,5,6,7,8,9,17,18,19,20,21,22,23,24),
                                       var.labels.reorder=c(12,13,2,1,11,7,8,9,10,3,4,5,6),
                                       color.panel=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7",
                                                     "#666666","#AD7700","#1C91D4","#D5C711","#A04700","#B14380","#4D4D4D","#FFBE2D"),
                                       scale="percent")+ggtitle("snRNA sample celltype proportions")+labs(x="Sample",y="% of cells")+
  mytheme+NoLegend()
ggsave('sample2_BarPlot_RNA.tiff', plot=sample2_BarPlot_RNA, height=80, width = 100, units="mm", dpi=300)

sample1_BarPlot_ATAC<-dittoBarPlot(intATAC, "cellsubtype",group.by="sample", 
                                        x.reorder=c(11,12,13,14,10,1,2,3,4,5,6,7,8,9,15,16,17,18,19,20,21,22), 
                                        ,var.labels.reorder=c(12,13,2,1,11,7,8,9,10,3,4,5,6),
                                        color.panel=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#CC79A7","#666666",
                                                      "#AD7700","#1C91D4","#A04700","#B14380","#4D4D4D","#FFBE2D"),
                                        scale="count")+ggtitle("snATAC sample cell numbers")+labs(x="Sample",y="Number of cells")+
  mytheme+NoLegend()
ggsave('sample1_BarPlot_ATAC.tiff', plot=sample1_BarPlot_ATAC, height=80, width = 100, units="mm", dpi=300)

sample2_BarPlot_ATAC<-dittoBarPlot(intATAC, "cellsubtype",group.by="sample", 
                                        x.reorder=c(11,12,13,14,10,1,2,3,4,5,6,7,8,9,15,16,17,18,19,20,21,22), 
                                        var.labels.reorder=c(12,13,2,1,11,7,8,9,10,3,4,5,6),
                                        color.panel=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#CC79A7","#666666",
                                                      "#AD7700","#1C91D4","#A04700","#B14380","#4D4D4D","#FFBE2D"),
                                        scale="percent")+ggtitle("snATAC sample celltype proportions")+labs(x="Sample",y="% of cells")+
  mytheme+NoLegend()
ggsave('sample2_BarPlot_ATAC.tiff', plot=sample2_BarPlot_ATAC, height=80, width = 100, units="mm", dpi=300)

p1a<-DimPlot_scCustom(itPFCn,group.by="cluster", colors_use =c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
  "#D55E00","#CC79A7","#666666","#AD7700","#1C91D4","#007756","#D5C711","#005685","#A04700","#B14380","#4D4D4D","#FFBE2D"))+
blank_theme+theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))+theme(title=element_text(size=7))+
guides(color = guide_legend(override.aes = list(size=2), ncol=1))

#ggsave("intRNA_batch_UMAPs.tiff", intRNA_batch_UMAPs, width=160, height=80, units='mm', dpi=300)
ggsave("reference_UMAP.pdf", p1a, width=100, height=80,units='mm')


p1b<-DimPlot_scCustom(intPFC,group.by="cluster", colors_use =c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2",
  "#D55E00","#CC79A7","#666666","#AD7700","#1C91D4","#007756","#D5C711","#005685","#A04700","#B14380","#4D4D4D","#FFBE2D"))+
blank_theme+theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))+theme(title=element_text(size=7))+NoLegend()
ggsave("reference_UMAP.tiff", p1b, width=80, height=80,units='mm')



get_cellsubtypes_legend<-DimPlot_scCustom(intRNA, group.by="Velm_cellsubtype",label=F,split.by="group",pt.size=0.0,raster=F,colors_use=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#CC79A7","#666666","#AD7700","#1C91D4","#A04700","#B14380","#4D4D4D","#FFBE2D"))+blank_theme+theme(legend.position="right", plot.title=element_blank(), legend.title = element_text(colour="black", size=7, face="bold"), legend.text=element_text(size=6), legend.key.size = unit(0.2, 'lines'))+
  guides(colour = guide_legend(override.aes = list(size = 4), ncol=1))
ggsave(get_cellsubtypes_legend, file="getlegend.pdf", height = 3.6, width = 7.2)

RNA_groups_cellsubtypes<-DimPlot_scCustom(intRNA, group.by="Velm_cellsubtype",label=F,split_seurat = TRUE, split.by="group",pt.size=0.001, colors_use=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#CC79A7","#666666","#AD7700","#1C91D4","#A04700","#B14380","#4D4D4D","#FFBE2D"))+blank_theme+theme(legend.position="right", plot.title=element_blank(), legend.title = element_text(colour="black", size=7, face="bold"), legend.text=element_text(size=5), legend.key.size = unit(0.2, 'lines'))+
guides(colour = guide_legend(override.aes = list(size = 4), ncol=2))
ggsave(RNA_groups_cellsubtypes, file="RNA_group_split_Dimplots.tiff", height = 3.6, width = 14.4)

ATAC_groups_cellsubtypes<-DimPlot_scCustom(intATAC, group.by="Velm_cellsubtype",label=F,split_seurat = TRUE, split.by="group",pt.size=0.001, colors_use = c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#CC79A7","#666666","#AD7700","#1C91D4","#A04700","#B14380","#4D4D4D","#FFBE2D"))+blank_theme+theme(legend.position="right", plot.title=element_blank(), legend.title = element_text(colour="black", size=7, face="bold"), legend.text=element_text(size=5), legend.key.size = unit(0.2, 'lines'))+
guides(colour = guide_legend(override.aes = list(size = 4), ncol=2))
ggsave(ATAC_groups_cellsubtypes, file="ATAC_group_split_Dimplots.tiff", height = 3.6, width = 14.4)


mytheme<-theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),plot.title = element_text(size=5, face="bold"),
               panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
               legend.title = element_text(colour="black", size=5, face="bold"),
               axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
               axis.title.x = element_text(colour="black", size=5, face="bold"), axis.title.y = element_text(colour="black", size=5, face="bold"))

###
diagnosis_Velm_BarPlot_RNA<-dittoBarPlot(intRNA, "Velm_cellsubtype",group.by="diagnosis",
x.reorder=c(4,3,1,2,5),var.labels.reorder=c(12,13,2,1,11,7,8,9,10,3,4,5,6),
color.panel=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#CC79A7",
"#666666","#AD7700","#1C91D4","#A04700","#B14380","#4D4D4D","#FFBE2D"),scale="count")+
ggtitle("snRNA diagnosis celltypes")+labs(x="Diagnosis",y="Number of cells")+
mytheme
#ggsave('figures/diagnosis_Velm_BarPlot_RNA.png', plot=diagnosis_Velm_BarPlot_RNA, width = 4.5, height = 4.5, dpi=300, units='in')
diagnosis_Velm_BarPlot_ATAC<-dittoBarPlot(intATAC, "Velm_cellsubtype",group.by="diagnosis",
x.reorder=c(4,3,1,2,5),var.labels.reorder=c(12,13,2,1,11,7,8,9,10,3,4,5,6),
color.panel=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#CC79A7",
"#666666","#AD7700","#1C91D4","#A04700","#B14380","#4D4D4D","#FFBE2D"), scale="count")+
ggtitle("snATAC diagnosis celltypes")+labs(x="Diagnosis",y="Number of cells")+
mytheme


diagnosis_Velm_BarPlot_RNA<-dittoBarPlot(intRNA, "Velm_cellsubtype",group.by="diagnosis",
x.reorder=c(4,3,1,2,5),var.labels.reorder=c(12,13,2,1,11,7,8,9,10,3,4,5,6),
color.panel=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#CC79A7",
"#666666","#AD7700","#1C91D4","#A04700","#B14380","#4D4D4D","#FFBE2D"),scale="percent")+
ggtitle("snRNA celltypes")+labs(x="Diagnosis",y="% of cells")+
mytheme+theme(legend.position="right", legend.text=element_text(size=5), legend.key.size = unit(0, 'lines'))+
  guides(color=guide_legend(ncol=1,title = "sample"))

diagnosis_Velm_BarPlot_ATAC<-dittoBarPlot(intATAC, "Velm_cellsubtype",group.by="diagnosis",
x.reorder=c(4,3,1,2,5),var.labels.reorder=c(12,13,2,1,11,7,8,9,10,3,4,5,6),
color.panel=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#CC79A7",
"#666666","#AD7700","#1C91D4","#A04700","#B14380","#4D4D4D","#FFBE2D"), scale="percent")+
ggtitle("snATAC celltypes")+labs(x="Diagnosis",y="% of cells")+
mytheme+theme(legend.position="right", legend.text=element_text(size=5), legend.key.size = unit(0, 'lines'))+
  guides(color=guide_legend(ncol=1,title = "sample"))

pdf("diagnosis_Velm_BarPlot_RNA.pdf",height = 2,width=2)
diagnosis_Velm_BarPlot_RNA
dev.off()
pdf("diagnosis_Velm_BarPlot_ATAC.pdf",height = 2,width=2)
diagnosis_Velm_BarPlot_ATAC
dev.off()

pseudobulk_Velm_BarPlot_merged<-dittoBarPlot(merger, "Velm_cellsubtype",group.by="pseudobulk",
scale="count",y.breaks=c(10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,110000,120000),
var.labels.reorder=c(12,13,2,1,11,7,8,9,10,3,4,5,6),
color.panel=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#CC79A7",
"#666666","#AD7700","#1C91D4","#A04700","#B14380","#4D4D4D","#FFBE2D"))+
ggtitle("Celltypes")+labs(x="Celltypes",y="Number of cells")+
mytheme+
mytheme+theme(legend.position="right", legend.text=element_text(size=5), legend.key.size = unit(0, 'lines'))+
guides(color=guide_legend(shape=15, ncol=1,title = "sample"))
pseudobulk_Velm_BarPlot_ATAC<-dittoBarPlot(intATAC, "Velm_cellsubtype",group.by="pseudobulk",scale="count",y.breaks=c(10000,20000,30000,40000,50000,60000,70000,80000,90000,100000),
var.labels.reorder=c(12,13,2,1,11,7,8,9,10,3,4,5,6),
color.panel=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#CC79A7",
"#666666","#AD7700","#1C91D4","#A04700","#B14380","#4D4D4D","#FFBE2D"))+
ggtitle("snATAC cells")+labs(x="Celltypes",y="Number of cells")+
mytheme+
mytheme+theme(legend.position="right", legend.text=element_text(size=5), legend.key.size = unit(0, 'lines'))+
guides(color=guide_legend(shape=15, ncol=1,title = "sample"))

pseudobulk_Velm_BarPlot_merged<-dittoBarPlot(merger, "Velm_cellsubtype",group.by="pseudobulk",
                                             scale="count",y.breaks=c(10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,110000,120000),
                                             x.reorder=c(2,1),var.labels.reorder=c(12,13,2,1,11,7,8,9,10,3,4,5,6),
                                             color.panel=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#CC79A7",
                                                           "#666666","#AD7700","#1C91D4","#A04700","#B14380","#4D4D4D","#FFBE2D"))+labs(x="Celltypes",y="Number of cells")+
  mytheme+  mytheme+theme(legend.position="right", legend.text=element_text(size=5), legend.key.size = unit(0, 'lines'),plot.title=element_blank())+
  guides(color=guide_legend(shape=15, ncol=1,title = "sample"))

pdf("pseudobulk_Velm_BarPlot_RNA.pdf",height = 2.2,width=2)
pseudobulk_Velm_BarPlot_RNA
dev.off()
pdf("pseudobulk_Velm_BarPlot_ATAC.pdf",height = 2.2,width=2)
pseudobulk_Velm_BarPlot_ATAC
dev.off()

p1<-DimPlot_scCustom(intRNA,group.by="diagnosis",pt.size=0.001, colors_use=group_colors)+blank_theme+theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))+theme(title=element_blank())+NoLegend()
p2<-DimPlot_scCustom(intATAC,group.by="diagnosis",pt.size=0.001, colors_use=group_colors)+blank_theme+theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"cm"))+theme(title=element_blank())
ggsave("diagnosis_plot.tiff",diag_plot,height=96, width=220,unit="mm")


RNA_markers <- c('MOBP','MYT1','AQP4','CD74','SLC17A7','GAD2')
pal <- viridis(n = 10, option = "D")
pList_RNAmarkers <- list()
for(i in 1:length(RNA_markers)){
  cur_marker<-RNA_markers[[i]]
  plot<-FeaturePlot_scCustom(intRNA,features=cur_marker)+blank_theme+theme(plot.title=element_blank())+NoLegend()+theme(plot.margin=unit(c(1,1.5,0,0),"cm"))
  pList_RNAmarkers[[i]] <- plot
}
snRNA_markers<-wrap_plots(pList_RNAmarkers)+plot_layout(ncol=3,nrow=2)
ggsave("snRNA_markers.tiff", snRNA_markers, width=240, height=160, units='mm', dpi=300)

ATAC_markers <- c('MOBP','MYT1','GFAP','CSF1R','NRGN','GAD2')
pal <- viridis(n = 10, option = "D")
pList_ATACmarkers <- list()
for(i in 1:length(ATAC_markers)){
  cur_marker<-ATAC_markers[[i]]
  plot<-FeaturePlot_scCustom(intATAC,features=cur_marker)+blank_theme+theme(plot.title=element_blank())+NoLegend()+theme(plot.margin=unit(c(1,1.5,0,0),"cm"))
  pList_ATACmarkers[[i]] <- plot
}
snATAC_markers<-wrap_plots(pList_ATACmarkers)+plot_layout(ncol=3,nrow=2)
ggsave("snATAC_markers.tiff", snATAC_markers, width=240, height=160, units='mm', dpi=300)


ATAC_markers <- c('MOBP','MYT1','GFAP','CSF1R','NRGN','GAD2')
pal <- viridis(n = 10, option = "D")
pList_ATACmarkers <- list()
for(i in 1:length(ATAC_markers)){
  cur_marker<-ATAC_markers[[i]]
  plot<-FeaturePlot_scCustom(intATAC,features=cur_marker)+blank_theme+theme(plot.title=element_blank())+NoLegend()+theme(plot.margin=unit(c(1,1.5,0,0),"cm"))
  pList_ATACmarkers[[i]] <- plot
}
snATAC_markers<-wrap_plots(pList_ATACmarkers)+plot_layout(ncol=3,nrow=2)
ggsave("snATAC_markers.tiff", snATAC_markers, width=240, height=160, units='mm', dpi=300)


motif_markers <- c('SOX13','SOX10','PAX4','SPIC','NEUROD2','ZBTB18')
pal <- viridis(n = 10, option = "D")
pList_motifmarkers <- list()
for(i in 1:length(motif_markers)){
  cur_marker<-motif_markers[[i]]
  plot<-FeaturePlot_scCustom(intATAC,features=cur_marker)+blank_theme+theme(plot.title=element_blank())+NoLegend()+theme(plot.margin=unit(c(1,1.5,0,0),"cm"))
  pList_motifmarkers[[i]] <- plot
}
snATAC_motif_markers<-wrap_plots(pList_motifmarkers)+plot_layout(ncol=3,nrow=2)
ggsave("snATAC_motifs.tiff", snATAC_motif_markers, width=240, height=160, units='mm', dpi=300)

pList_motifmarkers <- list()
for(i in 1:length(motif_markers)){
  cur_marker<-motif_markers[[i]]
  plot<-FeaturePlot_scCustom(intATAC,features=cur_marker)+blank_theme+theme(plot.title=element_blank())+theme(plot.margin=unit(c(1,1.5,0,0),"cm"))
  pList_motifmarkers[[i]] <- plot
}
snATAC_motif_markers<-wrap_plots(pList_motifmarkers)+plot_layout(ncol=3,nrow=2)
ggsave("snATAC_motifs.pdf", snATAC_motif_markers, width=240, height=160, units='mm')

# ALS risk genes
DefaultAssay(intRNA)<-"SCT"
table(intRNA$diagnosis_celltype)

Idents(intRNA)<-"diagnosis_celltype"

#compiled from known all known GWAS and risk genes
ALSgenes<-c("TARDBP","C9orf72","SOD1","FUS","NEK1","OPTN","CHCHD10","SQSTM1",
            "TBK1","KIF5A","SETX","UBQLN2","MATR3","VAPB","SIGMAR1","ANXA11",
            "TUBA4A","ALS2","GRN","PFN1","CHMP2B","TIA1","ANG","SPAST","FIG4",
            "SPG11","GLE1","CCNF","ATXN2","VCP")
ALSgenes_dotplot<-Clustered_DotPlot(intRNA,features=ALSgenes,x_lab_rotate=TRUE, k=7, seed=1234)

#risk genes from van Wheenen et al. Nature Genetics 2021
ALS_risk_genes<-c("MOBP","NEK1","TNIP1","GPX3","ERGIC1","HLA-A","PTPRN2",
  "C9orf72","KIF5A","TBK1","COG3","SCFD1","UNC13A","SOD1","CFAP410","SPATA2")
ALSgenes_dotplot<-Clustered_DotPlot(intRNA,features=ALS_risk_genes,x_lab_rotate=TRUE, k=6, seed=1234)


