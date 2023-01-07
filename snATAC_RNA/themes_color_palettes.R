library(Seurat)
library(Signac)
library(patchwork)
library(dittoSeq)
library(scCustomize)
library(ggplot2)
library(viridis)
set.seed(1234)
library(ggplot2)

###############
###color palettes
###############

diagnosis_colors<-c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2")
celltype_colors<-c("#D55E00","#CC79A7","#E69F00","#0072B2","#009F73","#F0E442")
cellsubtype_colors <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#CC79A7","#666666","#AD7700","#1C91D4","#A04700","#B14380","#4D4D4D","#FFBE2D")

correlation_celltype_colors<-c("#E69F00","#009F73","#F0E442","#0072B2","#D55E00","#CC79A7")
names(correlation_celltype_colors)<-celltypes

#sample palettes for each technology/pipeline
#24 samples
pal_snRNA<-c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7",
  "#666666","#AD7700","#1C91D4","#007756","#D5C711","#005685","#B14380","#4D4D4D",
  "#FFBE2D","#80C7EF","#00F6B3","#F4EB71","#06A5FF","#FF8320","#D99BBD","#8C8C8C","#FFCB57")
#22 samples
pal_snATAC<-c("#E69F00","#56B4E9","#009E73","#F0E442","#CC79A7","#666666",
  "#AD7700","#1C91D4","#D5C711","#005685","#A04700","#B14380","#4D4D4D","#FFBE2D",
  "#80C7EF","#00F6B3","#F4EB71","#06A5FF","#FF8320","#D99BBD","#8C8C8C","#FFCB57")

#for ArchR plots (order based on forced ArchR alphabetical)
pal_snATAC_ArchR<-c("#666666","#AD7700","#1C91D4","#D5C711","#005685","#A04700",
  "#B14380","#4D4D4D","#FFBE2D","#CC79A7","#E69F00","#56B4E9","#009E73","#F0E442",
  "#80C7EF","#00F6B3","#F4EB71","#06A5FF","#FF8320","#D99BBD","#8C8C8C","#FFCB57")

#for ArchR ridgeplots (needs reverse alphabetical)
pal_snATAC_ArchR2<-rev(pal_snATAC_ArchR)

#size 5 font theme
mytheme<-theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),plot.title = element_text(size=5, face="bold"),
panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
legend.title = element_text(colour="black", size=5, face="bold"),
axis.text.x = element_text(colour="black", size=5), axis.text.y = element_text(colour="black", size=5),
axis.title.x = element_text(colour="black", size=5, face="bold"), axis.title.y = element_text(colour="black", size=5, face="bold"))

#size 7 font theme
mytheme<-theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),plot.title = element_text(size=7, face="bold"),
panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
legend.title = element_text(colour="black", size=7, face="bold"),
axis.text.x = element_text(colour="black", size=7), axis.text.y = element_text(colour="black", size=7),
axis.title.x = element_text(colour="black", size=7, face="bold"), axis.title.y = element_text(colour="black", size=7, face="bold"))

blank_theme <- theme(
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
 # legend.position="none",
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  plot.background=element_blank()
)
