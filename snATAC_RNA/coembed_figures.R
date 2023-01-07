### coembed celltype vis
library(Seurat)
library(Signac)
library(ggplot2)
library(scCustomize)
set.seed(1234)

coembed <- readRDS("/mnt/WORKHORSE/Jan22_RNA_ATAC/Apr8_coembed.RDS")

coembed$Velm_cellsubtype <- factor(coembed$Velm_cellsubtype, levels = c("Oligodendrocytes","OPC","AST-PP","AST-FB","Microglia","Endothelial",
"L2/3","L4","L5/6","L5/6-CC","Neu-NRGN-II","IN-PV","IN-SST","IN-SV2C","IN-VIP"))
coembed$diagnosis<-factor(coembed$diagnosis, levels = c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD"))
coembed$celltype<-factor(coembed$celltype,levels=c("Oligodendrocytes","OPC","Astrocytes","Microglia","Endothelial","Excitatory","Inhibitory"))

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

group_colors<-c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2")
cellsubtype_colors <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#CC79A7","#666666","#AD7700","#1C91D4","#A04700","#B14380","#4D4D4D","#FFBE2D")
celltype_colors<-c("#D55E00","#CC79A7","#E69F00","#0072B2","#56B4E9","#009F73","#F0E442")

coembed_tech<-DimPlot_scCustom(coembed,group.by="technology", colors_use = safe_colorblind_palette,label=F)+
  theme(legend.position="right")+guides(color = guide_legend(override.aes = list(size=4), ncol=1))+ggtitle("coembed by technology")+labs(x="UMAP 1",y="UMAP 2")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        legend.title = element_text(colour="black", size=7, face="bold"),
        axis.text.x = element_text(colour="black", size=7), axis.text.y = element_text(colour="black", size=7),
        axis.title.x = element_text(colour="black", size=7, face="bold"), axis.title.y = element_text(colour="black", size=7, face="bold")) 
ggsave('coembed_tech.tiff', plot=coembed_tech, width = 100, height = 80, units='mm')

coembed_celltype<-DimPlot_scCustom(coembed,group.by="celltype", colors_use = celltype_colors,label=F)+
  theme(legend.position="right")+guides(color = guide_legend(override.aes = list(size=4), ncol=1))+ggtitle("coembed by major celltypes")+labs(x="UMAP 1",y="UMAP 2")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        legend.title = element_text(colour="black", size=7, face="bold"),
        axis.text.x = element_text(colour="black", size=7), axis.text.y = element_text(colour="black", size=7),
        axis.title.x = element_text(colour="black", size=7, face="bold"), axis.title.y = element_text(colour="black", size=7, face="bold")) 
ggsave('coembed_celltype.tiff', plot=coembed_celltype, width = 100, height = 80, units='mm')

coembed_Velm<-DimPlot_scCustom(coembed,group.by="Velm_cellsubtype", colors_use = cellsubtype_colors,label=F)+
  theme(legend.position="right")+guides(color = guide_legend(override.aes = list(size=4), ncol=1))+ggtitle("coembed by major celltypes")+labs(x="UMAP 1",y="UMAP 2")+
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
        legend.title = element_text(colour="black", size=7, face="bold"),
        axis.text.x = element_text(colour="black", size=7), axis.text.y = element_text(colour="black", size=7),
        axis.title.x = element_text(colour="black", size=7, face="bold"), axis.title.y = element_text(colour="black", size=7, face="bold")) 
ggsave('coembed_Velm.tiff', plot=coembed_Velm, width = 100, height = 80, units='mm')

options(scipen=5)

hist(coembed$prediction.score.max, main="Seurat Label Prediction Score Max",
     xlab="prediction score max",
     ylab="cell frequency",
     col="blue",
     border="gray",ylim=c(0,100000))

ggplot(coembed2@meta.data, aes(x=coembed2$prediction.score.max)) +
  geom_histogram( position="identity", alpha=0.5)+
  geom_vline(xintercept=0.5)
dev.off()


blank_theme <- theme(
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  legend.position="none",
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  plot.background=element_blank()
)
prediction_features <- c('prediction.score.Oligodendrocytes','prediction.score.OPC','prediction.score.AST.PP', 'prediction.score.AST.FB',
                      'prediction.score.Microglia','prediction.score.L2.3','prediction.score.L4','prediction.score.L5.6','prediction.score.L5.6.CC',
                      'prediction.score.IN.PV','prediction.score.IN.SST','prediction.score.IN.SV2C','prediction.score.IN.VIP')


prediction_features <- c(prediction_features, 'prediction.score.max')
plots <- FeaturePlot(coembed, features=prediction_features)
for(i in 1:length(plots)){
  plots[[i]] <- plots[[i]] + blank_theme
}
png(paste0('prediction_score_featureplot.tiff'), width=15, height=15, res=300, units='mm')
wrap_plots(plots, ncol=4)
dev.off()