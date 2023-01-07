intOPC$diagnosis<-intOPC$group
intOPC$diagnosis<-factor(intOPC$diagnosis,levels=diagnosis)
Idents(intOPC)<-'diagnosis'

NFE2L1_OPC<-VlnPlot_scCustom(intOPC,features="NFE2L1",group.by = "diagnosis",pt.size=0,colors_use = group_colors, assay="chromvar")+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("NFE2L1")+labs(x="",y="chromVAR deviations")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=5), axis.title.y = element_text(colour="black", size=5))
ggsave('NFE2L1_OPC.pdf', plot=NFE2L1_OPC, width = 100, height = 120, units='mm')

intEXC$diagnosis<-intEXC$group
intEXC$diagnosis<-factor(intEXC$diagnosis,levels=diagnosis)
Idents(intEXC)<-'diagnosis'

intINH$diagnosis<-intINH$group
intINH$diagnosis<-factor(intINH$diagnosis,levels=diagnosis)
Idents(intINH)<-'diagnosis'


RELA_EXC<-VlnPlot_scCustom(intEXC,features="RELA",group.by = "diagnosis",pt.size=0,colors_use = group_colors, assay="chromvar")+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("RELA")+labs(x="",y="chromVAR deviations")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=5), axis.title.y = element_text(colour="black", size=5))
ggsave('RELA_EXC.pdf', plot=RELA_EXC, width = 100, height = 120, units='mm')


NFYA_INH<-VlnPlot_scCustom(intINH,features="NFYA",group.by = "diagnosis",pt.size=0,colors_use = group_colors, assay="chromvar")+
  geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("NFYA")+labs(x="",y="chromVAR deviations")+
  theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
        axis.title.x = element_text(colour="black", size=5), axis.title.y = element_text(colour="black", size=5))
ggsave('NFYA_INH.pdf', plot=NFYA_INH, width = 100, height = 120, units='mm')


intAST$diagnosis<-intAST$group
intAST$diagnosis<-factor(intAST$diagnosis,levels=diagnosis)
Idents(intAST)<-'diagnosis'

SLC1A2_coverage<-CoveragePlot(intAST,region="SLC1A2",group.by="diagnosis")&
  scale_fill_manual(values=group_colors)&theme(axis.line = element_line(colour = "black"),
                                               panel.grid.major = element_blank(),plot.title = element_text(size=5, face="bold"),
                                               panel.grid.minor = element_blank(),panel.border = element_blank(),
                                               panel.background = element_blank(),
                                               legend.title = element_text(colour="black", size=5, face="bold"),
                                               axis.text.x = element_text(colour="black", size=5), 
                                               axis.text.y = element_text(colour="black", size=5),
                                               axis.title.x = element_text(colour="black", size=5, face="bold"), 
                                               axis.title.y = element_text(colour="black", size=5, face="bold"))&
  theme(panel.grid.major = element_line(colour="lightgray", size=0.25))

GFAP_coverage<-CoveragePlot(intAST,region="GFAP",group.by="diagnosis")&
  scale_fill_manual(values=group_colors)&theme(axis.line = element_line(colour = "black"),
                                               panel.grid.major = element_blank(),plot.title = element_text(size=5, face="bold"),
                                               panel.grid.minor = element_blank(),panel.border = element_blank(),
                                               panel.background = element_blank(),
                                               legend.title = element_text(colour="black", size=5, face="bold"),
                                               axis.text.x = element_text(colour="black", size=5), 
                                               axis.text.y = element_text(colour="black", size=5),
                                               axis.title.x = element_text(colour="black", size=5, face="bold"), 
                                               axis.title.y = element_text(colour="black", size=5, face="bold"))&
  theme(panel.grid.major = element_line(colour="lightgray", size=0.25))



atacMICRO<-readRDS(file="/mnt/WORKHORSE/Jan22_RNA_ATAC/objects_ATAC/Mar11_ATAC_microglia.RDS")
fragments<-CreateFragmentObject(path="/mnt/WORKHORSE/Jan22_RNA_ATAC/fragments/fragments.tsv.gz",cells=colnames(atacMICRO))
#intMICRO<-readRDS(file="~/projects/def-rogaeva/prime235/Oct21_ATAC/objects/Mar11_ATAC_microglia.RDS")
#fragments<-CreateFragmentObject(path="~/projects/def-rogaeva/prime235/Oct21_ATAC/fragments/fragments.tsv.gz",cells=colnames(intMICRO))
Fragments(atacMICRO)<-NULL
Fragments(atacMICRO)<-fragments


atacMICRO$diagnosis<-atacMICRO$group
atacMICRO$diagnosis<-factor(atacMICRO$diagnosis,levels=diagnosis)
Idents(atacMICRO)<-'diagnosis'

CoveragePlot(atacMICRO,region="CSF1R",group.by="diagnosis")&
  scale_fill_manual(values=group_colors)&theme(axis.line = element_line(colour = "black"),
                                               panel.grid.major = element_blank(),plot.title = element_text(size=5, face="bold"),
                                               panel.grid.minor = element_blank(),panel.border = element_blank(),
                                               panel.background = element_blank(),
                                               legend.title = element_text(colour="black", size=5, face="bold"),
                                               axis.text.x = element_text(colour="black", size=5), 
                                               axis.text.y = element_text(colour="black", size=5),
                                               axis.title.x = element_text(colour="black", size=5, face="bold"), 
                                               axis.title.y = element_text(colour="black", size=5, face="bold"))&
  theme(panel.grid.major = element_line(colour="lightgray", size=0.25))





raw_metrics <- Read_Metrics_10X(base_path = "/mnt/DATA1/Documents/single_cell/snRNAseq/cellranger5_snRNA_output/", default_10X = TRUE)


head(raw_metrics)
Seq_QC_Plot_Basic_Combined(metrics_dataframe = raw_metrics, plot_by = "sample_id")
Seq_QC_Plot_Alignment_Combined(metrics_dataframe = raw_metrics, plot_by = "sample_id")
Seq_QC_Plot_Basic_Combined(metrics_dataframe = raw_metrics, plot_by = "sample_id", significance = T)
polychrome_pal <- DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")
DimPlot(object = intRNA, cols = DiscretePalette_scCustomize(num_colors = 26, palette = "polychrome"))
polychrome_pal <- DiscretePalette_scCustomize(num_colors = 24, palette = "alphabet")
DimPlot(object = intRNA, cols = DiscretePalette_scCustomize(num_colors = 17, palette = "alphabet"))

c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD")
VlnPlot_scCustom(intAST, features="LHX9",colors_use = group_colors,group.by='group',pt.size=0)+geom_boxplot(width=0.2,fill="white")+ylab("chromVAR deviations")+NoLegend()+mytheme
VlnPlot_scCustom(intAST, features="PAX4",colors_use = group_colors,group.by='group',pt.size=0)+geom_boxplot(width=0.2,fill="white")+ylab("chromVAR deviations")+NoLegend()+mytheme
VlnPlot_scCustom(intMICRO, features="SPI1",colors_use = group_colors,group.by='group',pt.size=0)+geom_boxplot(width=0.2,fill="white")+ylab("chromVAR deviations")+NoLegend()+mytheme
VlnPlot_scCustom(intMICRO, features="SPIC",colors_use = group_colors,group.by='group',pt.size=0)+geom_boxplot(width=0.2,fill="white")+ylab("chromVAR deviations")+NoLegend()+mytheme
VlnPlot_scCustom(intEXC, features="NEUROD2",colors_use = group_colors,group.by='group',pt.size=0)+geom_boxplot(width=0.2,fill="white")+ylab("chromVAR deviations")+NoLegend()+mytheme
VlnPlot_scCustom(intINH, features="ZBTB18",colors_use = group_colors,group.by='group',pt.size=0)+geom_boxplot(width=0.2,fill="white")+ylab("chromVAR deviations")+NoLegend()+mytheme

VlnPlot_scCustom(intMICRO, features="CTCF",colors_use = group_colors,group.by='group',pt.size=0)+geom_boxplot(width=0.2,fill="white")+ylab("chromVAR deviations")+NoLegend()+mytheme

C9orf72 coverageplots


sample1_Velm_BarPlot_RNA<-dittoBarPlot(intRNA, "Velm_cellsubtype",group.by="sample", 
                                       x.reorder=c(11,12,13,14,15,16,10,1,2,3,4,5,6,7,8,9,17,18,19,20,21,22,23,24),
                                       var.labels.reorder=c(12,13,2,1,11,7,8,9,10,3,4,5,6),
                                       color.panel=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7",
                                                     "#666666","#AD7700","#1C91D4","#D5C711","#A04700","#B14380","#4D4D4D","#FFBE2D"),
                                       scale="count")+ggtitle("snRNA sample cell numbers")+labs(x="Sample",y="Number of cells")+
  mytheme+NoLegend()
#ggsave('sample1_Velm_BarPlot_RNA.pdf', plot=sample1_Velm_BarPlot_RNA, width = 6, height = , units='mm')


sample2_Velm_BarPlot_RNA<-dittoBarPlot(intRNA, "Velm_cellsubtype",group.by="sample", 
                                       x.reorder=c(11,12,13,14,15,16,10,1,2,3,4,5,6,7,8,9,17,18,19,20,21,22,23,24),
                                       var.labels.reorder=c(12,13,2,1,11,7,8,9,10,3,4,5,6),
                                       color.panel=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7",
                                                     "#666666","#AD7700","#1C91D4","#D5C711","#A04700","#B14380","#4D4D4D","#FFBE2D"),
                                       scale="percent")+ggtitle("snRNA sample celltype proportions")+labs(x="Sample",y="% of cells")+
  mytheme+NoLegend()
#ggsave('sample2_Velm_BarPlot_RNA.tiff', plot=sample2_Velm_BarPlot_RNA, height=80, width = 100, units="mm", dpi=300)

sample1_Velm_BarPlot_ATAC<-dittoBarPlot(intATAC, "Velm_cellsubtype",group.by="sample", 
                                        x.reorder=c(11,12,13,14,10,1,2,3,4,5,6,7,8,9,15,16,17,18,19,20,21,22), 
                                        ,var.labels.reorder=c(12,13,2,1,11,7,8,9,10,3,4,5,6),
                                        color.panel=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#CC79A7","#666666",
                                                      "#AD7700","#1C91D4","#A04700","#B14380","#4D4D4D","#FFBE2D"),
                                        scale="count")+ggtitle("snATAC sample cell numbers")+labs(x="Sample",y="Number of cells")+
  mytheme+NoLegend()
#ggsave('sample1_Velm_BarPlot_ATAC.tiff', plot=sample1_Velm_BarPlot_ATAC, height=80, width = 100, units="mm", dpi=300)

sample2_Velm_BarPlot_ATAC<-dittoBarPlot(intATAC, "Velm_cellsubtype",group.by="sample", 
                                        x.reorder=c(11,12,13,14,10,1,2,3,4,5,6,7,8,9,15,16,17,18,19,20,21,22), 
                                        var.labels.reorder=c(12,13,2,1,11,7,8,9,10,3,4,5,6),
                                        color.panel=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#CC79A7","#666666",
                                                      "#AD7700","#1C91D4","#A04700","#B14380","#4D4D4D","#FFBE2D"),
                                        scale="percent")+ggtitle("snATAC sample celltype proportions")+labs(x="Sample",y="% of cells")+
  mytheme+NoLegend()
#ggsave('sample2_Velm_BarPlot_ATAC.tiff', plot=sample2_Velm_BarPlot_ATAC, height=80, width = 100, units="mm", dpi=300)



CoveragePlot(intATAC,group.by="diagnosis",region="chr9-27325209-27642852",assay='peaks')&scale_fill_manual(values=group_colors)+mytheme

blank_theme <- theme(
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  legend.position="right",
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  plot.background=element_blank()
)
pal<- viridis(n=10, option="D")

sex_DimPlot_RNA<-DimPlot_scCustom(intRNA, group.by="sex",colors_use=safe_colorblind_palette)+blank_theme+theme(text=element_text(size=12,family="Arial"))+ggtitle("snRNA by Sex")
sex_DimPlot_ATAC<-DimPlot_scCustom(intATAC, group.by="sex",colors_use=safe_colorblind_palette)+blank_theme+theme(text=element_text(size=12,family="Arial"))+ggtitle("snATAC by Sex")
ggsave("sex_DimPlot_RNA.tiff", sex_DimPlot_RNA, height=80, width = 100, units="mm", dpi=300)
ggsave("sex_DimPlot_ATAC.tiff", sex_DimPlot_ATAC, height=80, width = 100, units="mm", dpi=300)


chemistry_DimPlot_RNA<-DimPlot_scCustom(intRNA,group.by="chemistry")+theme(legend.position="right")+blank_theme+theme(text=element_text(size=12,family="Arial"))+
  ggtitle("snRNA by 10X chemistry")+labs(x="UMAP 1",y="UMAP 2")+
  blank_theme
ggsave("chemistry_DimPlot_RNA.tiff", chemistry_DimPlot_RNA, height=80, width = 100, units="mm", dpi=300)


nFeature_plot_RNA<-FeaturePlot_scCustom(intRNA,features="nFeature_RNA", colors_use = pal)+blank_theme+theme(text=element_text(size=12,family="Arial"))+
  scale_fill_continuous(limits=c())+ggtitle("snRNA number of genes")+labs(x="UMAP 1",y="UMAP2")+
  blank_theme
ggsave('nFeature_FeaturePlot_RNA.tiff', plot=nFeature_plot_RNA, height=80, width = 100, units="mm", dpi=300)

nCount_plot_RNA<-FeaturePlot_scCustom(intRNA,features="nCount_RNA", na_cutoff=1e-09,colors_use = pal)+blank_theme+theme(text=element_text(size=12,family="Arial"))+
  ggtitle("snRNA UMI counts (RNA)")+labs(x="UMAP 1",y="UMAP 2")
ggsave('nCount_FeaturePlot_RNA.tiff', plot=nCount_plot_RNA, height=80, width = 100, units="mm", dpi=300)

#Feature Plots SCT
nFeature_plot_SCT<-FeaturePlot_scCustom(intRNA,features="nFeature_SCT", colors_use = pal)+
  ggtitle("snRNA Features (SCT)")+blank_theme+theme(text=element_text(size=12,family="Arial"))
ggsave('nFeature_FeaturePlot_SCT.tiff', plot=nFeature_plot_SCT, height=80, width = 100, units="mm", dpi=300)

nCount_plot_SCT<-FeaturePlot_scCustom(intRNA,features="nCount_SCT", na_cutoff=1e-09,colors_use = pal)+
  ggtitle("snRNA UMI counts (SCT)")+blank_theme+theme(text=element_text(size=12,family="Arial"))
ggsave('nCount_FeaturePlot_SCT.tiff', plot=nCount_plot_SCT, height=80, width = 100, units="mm", dpi=300)

mito_plot_RNA<-FeaturePlot_scCustom(intRNA,features="percent.mt", colors_use = pal)+blank_theme+theme(text=element_text(size=12,family="Arial"))+
  ggtitle("snRNA %mitocondrial reads")+labs(x="UMAP 1",y="UMAP 2")+
  blank_theme
ggsave('mito_FeaturePlot_RNA.tiff', plot=mito_plot_RNA, height=80, width = 100, units="mm", dpi=300)

rpl_plot_RNA<-FeaturePlot_scCustom(intRNA,features="percent.rpl", colors_use = pal)+blank_theme+theme(text=element_text(size=12,family="Arial"))+
  ggtitle("snRNA %RPL reads")+labs(x="UMAP 1",y="UMAP 2")
ggsave('RPL_FeaturePlot_RNA.tiff', plot=rpl_plot_RNA, height=80, width = 100, units="mm", dpi=300)

rps_plot_RNA<-FeaturePlot_scCustom(intRNA,features="percent.rps", colors_use = pal)+blank_theme+theme(text=element_text(size=12,family="Arial"))+
  ggtitle("snRNA %RPS reads")+labs(x="UMAP 1",y="UMAP 2")
ggsave('RPS_FeaturePlot_RNA.tiff', plot=rps_plot_RNA, height=80, width = 100, units="mm", dpi=300)

Sscore_plot_RNA<-FeaturePlot_scCustom(intRNA,features="S.Score", colors_use = pal)+blank_theme+theme(text=element_text(size=12,family="Arial"))+
  ggtitle("snRNA %S phase reads")+labs(x="UMAP 1",y="UMAP 2")
ggsave('Sscore_FeaturePlot_RNA.tiff', plot=Sscore_plot_RNA, height=80, width = 100, units="mm", dpi=300)

G2Mscore_plot_RNA<-FeaturePlot_scCustom(intRNA,features="G2M.Score", colors_use = pal)+blank_theme+theme(text=element_text(size=12,family="Arial"))+
  ggtitle("snRNA %G2M phase reads")+labs(x="UMAP 1",y="UMAP 2")
ggsave('G2Mscore_FeaturePlot_RNA.tiff', plot=G2Mscore_plot_RNA, height=80, width = 100, units="mm", dpi=300)

CCphase_plot_RNA<-DimPlot_scCustom(intRNA, group.by="Phase")+blank_theme+theme(text=element_text(size=12,family="Arial"))+
  ggtitle("snRNA Cell Cycle Phase")+labs(x="UMAP 1",y="UMAP 2")
ggsave('CCphase_DimPlot_RNA.tiff', plot=CCphase_plot_RNA, height=80, width = 100, units="mm", dpi=300)


#sALSnoFTLD vs control

prediction_features <- c('prediction.score.Oligodendrocytes','prediction.score.OPC','prediction.score.AST.PP', 'prediction.score.AST.FB',
                         'prediction.score.Microglia','prediction.score.L2.3','prediction.score.L4','prediction.score.L5.6','prediction.score.L5.6.CC',
                         'prediction.score.IN.PV','prediction.score.IN.SST','prediction.score.IN.SV2C','prediction.score.IN.VIP','prediction.score.max')
prediction_features 

pList <- list ()
for(i in 1:length(prediction_features)){
  cur_prediction <- prediction_features[[i]]
  plot <- FeaturePlot(coembed,features=cur_prediction,raster=F) +my_theme+theme(text=element_text(size=12,family="Arial"))
  print(plot)
  pList[[i]] <- plot
}
coembed_all<-wrap_plots(pList,ncol=4)
plot <- FeaturePlot(coembed,features=prediction.score.max,raster=F) +my_theme+theme(text=element_text(size=12,family="Arial"))

png(paste0('prediction_score_featureplot.tiff'), width=30, height=30, res=300, units='cm')
wrap_plots(plots, ncol=4)
dev.off()

dev.off()
BiocManager::install("ReactomeGSA")
library(Reactome.GSA)
library(ReactomeGSA)
intMICRO
gsva_result<-analyse_sc_clusters(intMICRO, verbose=T)
gsava_result
gsva_result
pathway_expression<-pathways(gsva_result)

colnaems(pathways_epxression<-gsub("\\.Seurat","",colnames(pathway_expression)))
colnames(pathways_epxression<-gsub("\\.Seurat","",colnames(pathway_expression)))
pathway_expression[1:3,]
max_difference<-do.call(rbind,apply(pathway_expression,1,function(row){
  values<-as.numeric(row[2:length(row)])
  return(data.frame(name=row[1],min=min(values),max=max(values)))
}))
max_difference$diff<-max_difference$max-max_difference$min
max_difference<-max_difference[order(max_difference$diff,decreasing = T), ]
head(max_difference)
plot_gsva_pathway(gsva_result,pathays_id=rownames(max_difference[1])
                  plot_gsva_pathway(gsva_result,pathays_id=rownames(max_difference[1]))
                  plot_gsva_pathway(gsva_result,pathway_id=rownames(max_difference[1]))
                  plot_gsva_pathway(gsva_result,pathway_id=rownames(max_difference[2]))
                  plot_gsva_pathway(gsva_result,pathway_id=rownames(max_difference[3]))
                  max_difference
                  plot_gsva_heatmap(gsva_result, max_pathways = 15, margins = c(6,20))
                  dev.off()
                  plot_gsva_heatmap(gsva_result, max_pathways = 15, margins = c(6,20))
                  plot_gsva_pca(gsva_result)
                  
                  
                  blank_theme <- theme(
                    #text=element_blank(),
                    plot.title = element_blank(),
                    axis.line=element_blank(),
                    axis.text.x=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks=element_blank(),
                    axis.title.x=element_blank(),
                    axis.title.y=element_blank(),
                    legend.position="right",
                    panel.background=element_blank(),
                    panel.border=element_blank(),
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    plot.background=element_blank()
                  )
                  pal<- viridis(n=10, option="D")
                  
                  diagnosis_C9_Coverage<-CoveragePlot(intATAC,region="C9orf72",group.by="diagnosis")&scale_fill_manual(values=group_colors)&theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),plot.title = element_text(size=5, face="bold"),
                                                                                                                                                  panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
                                                                                                                                                  legend.title = element_text(colour="black", size=5, face="bold"),
                                                                                                                                                  axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
                                                                                                                                                  axis.title.x = element_text(colour="black", size=5, face="bold"), axis.title.y = element_text(colour="black", size=5, face="bold"))&
                    theme(panel.grid.major = element_line(colour="lightgray", size=0.25))
                  
                  DefaultAssay(intATAC)<-'peaks'
                  diagnosis_C9_Coverage<-CoveragePlot(intATAC,region="C9orf72",group.by="celltype")&scale_fill_manual(values=celltype_colors)&theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),plot.title = element_text(size=5, face="bold"),
                                                                                                                                                    panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
                                                                                                                                                    legend.title = element_text(colour="black", size=5, face="bold"),
                                                                                                                                                    axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
                                                                                                                                                    axis.title.x = element_text(colour="black", size=5, face="bold"), axis.title.y = element_text(colour="black", size=5, face="bold"))&
                    theme(panel.grid.major = element_line(colour="lightgray", size=0.25))
                  Idents(intATAC)<-'pseudobulk'
                  table(Idents(intATAC))
                  
                  DefultAssay(intRNA)<-'SCT'
                  DefaultAssay(intATAC)<-"RNA"
                  
                  
                  C9_FeaturePlot_RNA<-FeaturePlot_scCustom(intRNA,features="C9orf72",colors_use=viridis_inferno_dark_high)+blank_theme+NoLegend()+theme(plot.margin=unit(c(-0.3,-0.3,-0.3,-0.3),"cm"))
                  C9_FeaturePlot_ATAC<-FeaturePlot_scCustom(intATAC,features="C9orf72",colors_use=viridis_inferno_dark_high)+blank_theme+NoLegend()+theme(plot.margin=unit(c(-0.3,-0.3,-0.3,-0.3),"cm"))
                  ggsave("C9_FeaturePlot_RNA.tiff", C9_FeaturePlot_RNA, height=80, width = 100, units="mm", dpi=300)
                  ggsave("C9_FeaturePlot_ATAC.tiff", C9_FeaturePlot_ATAC, height=80, width = 100, units="mm", dpi=300)
                  
                  
                  
                  blank_theme <- theme(
                    #text=element_blank(),
                    plot.title = element_blank(),
                    axis.line=element_blank(),
                    axis.text.x=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks=element_blank(),
                    axis.title.x=element_blank(),
                    axis.title.y=element_blank(),
                    legend.position="right",
                    panel.background=element_blank(),
                    panel.border=element_blank(),
                    panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),
                    plot.background=element_blank()
                  )
                  
                  C9_FeaturePlot_RNA<-FeaturePlot_scCustom(intRNA,features="C9orf72",colors_use=viridis_inferno_dark_high)+blank_theme
                  C9_FeaturePlot_ATAC<-FeaturePlot_scCustom(intATAC,features="C9orf72",colors_use=viridis_inferno_dark_high)+blank_theme
                  ggsave("C9_FeaturePlot_RNA.pdf", C9_FeaturePlot_RNA, height=80, width = 100, units="mm")
                  ggsave("C9_FeaturePlot_ATAC.pdf", C9_FeaturePlot_ATAC, height=80, width = 100, units="mm")
                  
                  
                  #read and add fragment path
                  intATAC<-readRDS(file="~/projects/def-rogaeva/prime235/Oct21_ATAC/objects/Mar11_ATAC_final.RDS")
                  fragments<-CreateFragmentObject(path="~/projects/def-rogaeva/prime235/Oct21_ATAC/fragments/fragments.tsv.gz",cells=colnames(intATAC))
                  Fragments(intATAC)<-NULL
                  Fragments(intATAC)<-fragments

                   intATAC<-readRDS(file="/mnt/WORKHORSE/Jan22_RNA_ATAC/objects_ATAC/Mar11_ATAC_final.RDS")
                  fragments<-CreateFragmentObject(path="/mnt/WORKHORSE/Jan22_RNA_ATAC/fragments/fragments.tsv.gz",cells=colnames(intATAC))
                  Fragments(intATAC)<-NULL
                  Fragments(intATAC)<-fragments

                  intOLIGO$diagnosis<-intOLIGO$group
                  intOLIGO$diagnosis<-factor(intOLIGO$diagnosis,levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))
                  Idents(intOLIGO)<-"diagnosis"
                  
                  #intEXC<-readRDS(file="~/projects/def-rogaeva/prime235/Oct21_ATAC/objects/Mar11_ATAC_excitatory.RDS")
                  #fragments<-CreateFragmentObject(path="~/projects/def-rogaeva/prime235/Oct21_ATAC/fragments/fragments.tsv.gz",cells=colnames(intEXC))
                  Fragments(intOLIGO)<-NULL
                  Fragments(intOLIGO)<-fragments
                  intOLIGO$diagnosis<-intOLIGO$group
                  intOLIGO$diagnosis<-factor(intOLIGO$diagnosis,levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))
                  Idents(intOLIGO)<-"diagnosis"
                  CoveragePlot(intOLIGO,region="C9orf72",group.by="diagnosis")&scale_fill_manual(values=group_colors)&theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),plot.title = element_text(size=5, face="bold"),
                                                                                                                            panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
                                                                                                                            legend.title = element_text(colour="black", size=5, face="bold"),
                                                                                                                            axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
                                                                                                                            axis.title.x = element_text(colour="black", size=5, face="bold"), axis.title.y = element_text(colour="black", size=5, face="bold"))&
                    theme(panel.grid.minor = element_line(colour="lightgray", size=0.25))
                  
                  intAST$diagnosis<-intAST$group
                  intAST$diagnosis<-factor(intAST$diagnosis,levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))
                  Idents(intAST)<-"diagnosis"
                  
                  
                  CEBPG<-VlnPlot_scCustom(intMICRO,features="CEBPG",group.by = "group",pt.size=0,colors_use = group_colors)+
                    geom_boxplot(width=0.2,fill="white")+NoLegend()+ggtitle("CEBPg")+labs(x="",y="chromVAR deviations")+
                    theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
                          axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
                          axis.title.x = element_text(colour="black", size=6), axis.title.y = element_text(colour="black", size=6))
                  ggsave('CEBPG.pdf', plot=CEBPG, width = 3, height = 5, units='in')
                  
                  
                  intOPC$diagnosis<-intOPC$group
                  intOPC$diagnosis<-factor(intOPC$diagnosis,levels=c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))
                  Idents(intOPC)<-"diagnosis"
                  CoveragePlot(intOPC,region="C9orf72",group.by="diagnosis")&scale_fill_manual(values=group_colors)&theme(axis.line = element_line(colour = "black"),
                                                                                                                          panel.grid.major = element_blank(),plot.title = element_text(size=5, face="bold"),panel.grid.minor = element_blank(),
                                                                                                                          panel.border = element_blank(),panel.background = element_blank(),legend.title = element_text(colour="black", size=5, face="bold"),
                                                                                                                          axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),axis.title.x = element_text(colour="black", size=5, face="bold"), 
                                                                                                                          axis.title.y = element_text(colour="black", size=5, face="bold"))&theme(panel.grid.minor = element_line(colour="lightgray", size=0.25))
                  
                  library(readr)
                  EDLEE_up <- read_csv("ref_data/EDLEE_2019_UP.csv")
                  EDLEE_up
                  EDLEE_down <- read_csv("ref_data/EDLEE_2019_DOWN.csv")
                  EDLEE_down
                  MA_2022 <- read_csv("ref_data/MA_2022.csv")
                  MA_2022<-deframe(MA_2022)
                  MA_2022
                  
                  Idents(intL23)<-'diagnosis'
                  Idents(intL4)<-'diagnosis'
                  Idents(intPV)<-'diagnosis'
                  Idents(intSST)<-'diagnosis'
                  Idents(intSV2C)<-'diagnosis'
                  Idents(intVIP)<-'diagnosis'
                  
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
                  
                  ma22<
                    ma22<-as.data.frame(MA_2022)
                  ma22
                  ma22_df<-data.frame()
                  row.names(ma22_df)<-MA_2022
                  
                  idx = which(x %in% y) # Positions of the values of y in x
                  x = x[-idx] # Remove those values using their position and "-" operator
                  
                  library(dittoSeq)
                  ex_heatmap_C9vsCTRL<-dittoHeatmap(intRNA.L23.avg,genes=list(tibble(MA_2022)),main = "Excitatory DE MA2022 in C9ALSFTLD",
                                                    annot.by="diagnosis",annot.colors = group_colors, show_colnames = F, show_rownames = F)
                  ggsave('figures/ex_heatmap_C9vsCTRL.png', plot=ex_heatmap_C9vsCTRL, width = 5, height = 6, units='in')
                  
                  DefaultAssay(intINH)<-"SCT"
                  
                  intL23<-AddModuleScore(intL23, features=MA_2022, name = "MA2022",ctrl=5)
                  intINH<-AddModuleScore(intINH, features=ribo_inhibit, name = "ribogenes",ctrl=5)
                  
                  
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
                  
                  
                  inhSST_C9_heatmap<-dittoHeatmap(intRNA.INSST.avg,genes=DE_TDP_INSST,group.by="diagnosis", annot.by="diagnosis",highlight.features=c("overlap2_INSST"),complex=F)
                  
                  
                  
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
                  
                  
                  DE_OLIGO<-cellsubtype.C9ALSFTLD.cutoff$gene[cellsubtype.C9ALSFTLD.cutoff$Velm_cellsubtype == "Oligodendrocytes"]
                  DE_OPC<-cellsubtype.C9ALSFTLD.cutoff$gene[cellsubtype.C9ALSFTLD.cutoff$Velm_cellsubtype == "OPC"]
                  DE_ASTPP<-cellsubtype.C9ALSFTLD.cutoff$gene[cellsubtype.C9ALSFTLD.cutoff$Velm_cellsubtype == "AST-PP"]
                  DE_micro<-cellsubtype.C9ALSFTLD.cutoff$gene[cellsubtype.C9ALSFTLD.cutoff$Velm_cellsubtype == "Microglia"]
                  
                  OLIGO_C9_heatmap<-dittoHeatmap(intRNA.OLIGO.avg,genes=DE_OLIGO[1:50],main = "Oligo DE in C9ALSFTLD",
                                                 annot.by="diagnosis",annot.colors = group_colors, show_colnames = F, show_rownames = T,border_color = "lightgray")
                  #ggsave('figures/micro_heatmap_C9vsCTRL.pdf', plot=micro_C9_heatmap, width = 5, height = 6, units='in')
                  
                  OPC_C9_heatmap<-dittoHeatmap(intRNA.OPC.avg,genes=DE_OPC[1:50],main = "OPC DE in C9ALSFTLD",
                                               annot.by="diagnosis",annot.colors = group_colors, show_colnames = F, show_rownames = T,border_color = "lightgray")
                  #ggsave('figures/micro_heatmap_C9vsCTRL.pdf', plot=micro_C9_heatmap, width = 5, height = 6, units='in')
                  
                  ASTPP_C9_heatmap<-dittoHeatmap(intRNA.micro.avg,genes=DE_ASTPP[1:50],main = "ASTPP DE in C9ALSFTLD",
                                                 annot.by="diagnosis",annot.colors = group_colors, show_colnames = F, show_rownames = T,border_color = "lightgray")
                  #ggsave('figures/micro_heatmap_C9vsCTRL.pdf', plot=micro_C9_heatmap, width = 5, height = 6, units='in')
                  
                  micro_C9_heatmap<-dittoHeatmap(intRNA.micro.avg,genes=DE_micro[1:50],main = "micro DE in C9ALSFTLD",
                                                 annot.by="diagnosis",annot.colors = group_colors, show_colnames = F, show_rownames = T,border_color = "lightgray")
                  #ggsave('figures/micro_heatmap_C9vsCTRL.pdf', plot=micro_C9_heatmap, width = 5, height = 6, units='in')
                  
                  
                  ##Protoplasmic astrocytes
                  DEe_C9_AST<-DEenrichRPlot(
                    intASTPP,
                    ident.1 = "C9ALSFTLD",
                    ident.2 = "control",
                    balanced = TRUE,
                    logfc.threshold = 0.1,
                    assay = "RNA",
                    test.use = "LR",
                    latent.vars = c("sex","nCount_RNA"),
                    p.val.cutoff = 0.01,
                    max.genes=Inf,
                    enrich.database = "GO_Biological_Process_2021",
                    num.pathway = 25,
                    return.gene.list=TRUE
                  )
                  
                  DEe_C9_AST$pos$term<-gsub('.{13}$','',DEe_C9_AST$pos$term)
                  DEe_C9_AST$neg$term<-gsub('.{13}$','',DEe_C9_AST$neg$term)
                  
                  p1 <- ggplot(data = DEe_C9_AST$pos, aes_string(x=reorder(DEe_C9_AST$pos$term,DEe_C9_AST$pos$log10pval), y="log10pval")) +
                    geom_bar(stat = "identity", fill = "indianred3") +
                    coord_flip() + xlab("Pathway") +
                    ylab("-log10(pval)") +
                    ggtitle(paste(ident.1, "- enriched", enrich.database, "terms")) +
                    theme_classic() +
                    geom_text(aes_string(label = "term", y = 0),
                              size = 3.5,
                              color = "black",
                              position = position_dodge(1),
                              hjust = 0)+
                    theme(axis.title.y= element_blank(),
                          axis.text.y = element_blank(),
                          axis.ticks.y = element_blank())
                  p1
                  
                  p2 <- ggplot(data = DEe_C9_AST$neg, aes_string(x=reorder(DEe_C9_AST$neg$term,DEe_C9_AST$neg$log10pval), y = "log10pval")) +
                    geom_bar(stat = "identity", fill = "dodgerblue") +
                    coord_flip() + xlab("Pathway") +
                    ylab("-log10(pval)") +
                    ggtitle(paste(ident.1, "- depleted", enrich.database, "terms")) +
                    theme_classic() +
                    geom_text(aes_string(label = "term", y = 0),
                              size = 3.5,
                              color = "black",
                              position = position_dodge(1),
                              hjust = 0)+
                    theme(axis.title.y= element_blank(),
                          axis.text.y = element_blank(),
                          axis.ticks.y = element_blank())
                  p <- p1 + p2
                  p
                  
                  
                  ##OPC
                  DEe_C9_OPC<-DEenrichRPlot(
                    intOPC,
                    ident.1 = "C9ALSFTLD",
                    ident.2 = "control",
                    balanced = TRUE,
                    logfc.threshold = 0.1,
                    assay = "RNA",
                    test.use = "LR",
                    latent.vars = c("sex","nCount_RNA"),
                    p.val.cutoff = 0.01,
                    max.genes=Inf,
                    enrich.database = "GO_Biological_Process_2021",
                    num.pathway = 25,
                    return.gene.list=TRUE
                  )
                  
                  DEe_C9_OPC$pos$term<-gsub('.{13}$','',DEe_C9_OPC$pos$term)
                  DEe_C9_OPC$neg$term<-gsub('.{13}$','',DEe_C9_OPC$neg$term)
                  
                  p3 <- ggplot(data = DEe_C9_OPC$pos, aes_string(x=reorder(DEe_C9_OPC$pos$term,DEe_C9_OPC$pos$log10pval), y="log10pval")) +
                    geom_bar(stat = "identity", fill = "indianred3") +
                    coord_flip() + xlab("Pathway") +
                    ylab("-log10(pval)") +
                    ggtitle(paste(ident.1, "- enriched", enrich.database, "terms")) +
                    theme_classic() +
                    geom_text(aes_string(label = "term", y = 0),
                              size = 3.5,
                              color = "black",
                              position = position_dodge(1),
                              hjust = 0)+
                    theme(axis.title.y= element_blank(),
                          axis.text.y = element_blank(),
                          axis.ticks.y = element_blank())
                  
                  
                  p4 <- ggplot(data = DEe_C9_OPC$neg, aes_string(x=reorder(DEe_C9_OPC$neg$term,DEe_C9_OPC$neg$log10pval), y = "log10pval")) +
                    geom_bar(stat = "identity", fill = "dodgerblue") +
                    coord_flip() + xlab("Pathway") +
                    ylab("-log10(pval)") +
                    ggtitle(paste(ident.1, "- depleted", enrich.database, "terms")) +
                    theme_classic() +
                    geom_text(aes_string(label = "term", y = 0),
                              size = 3.5,
                              color = "black",
                              position = position_dodge(1),
                              hjust = 0)+
                    theme(axis.title.y= element_blank(),
                          axis.text.y = element_blank(),
                          axis.ticks.y = element_blank())
                  p5 <- p3 + p4
                  p5