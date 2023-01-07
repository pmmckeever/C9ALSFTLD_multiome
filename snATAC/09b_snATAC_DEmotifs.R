library(Seurat)
library(Signac)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dittoSeq)
library(BiocParallel)
register(MulticoreParam(48)) 
set.seed(1234)


diagnosis_colors<-c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2")
diagnosis<-c("control","C9noALSnoFTLD","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD")

#read and add fragment path
intATAC <- readRDS("objects/intATAC_final.RDS")
fragments<-CreateFragmentObject(path="fragments/fragments.tsv.gz",cells=colnames(intATAC))
Fragments(intATAC)<-NULL
Fragments(intATAC)<-fragments

#read and add fragment path
intOLIGO <- readRDS("objects/ATAC_oligo.RDS")
fragments<-CreateFragmentObject(path="fragments/fragments.tsv.gz",cells=colnames(intOLIGO))
Fragments(intOLIGO)<-NULL
Fragments(intOLIGO)<-fragments

#read and add fragment path
intOPC<-readRDS(file="objects/ATAC_OPC.RDS")
fragments<-CreateFragmentObject(path="fragments/fragments.tsv.gz",cells=colnames(intOPC))
Fragments(intOPC)<-NULL
Fragments(intOPC)<-fragments

#read and add fragment path
intAST<-readRDS(file="objects/ATAC_astrocytes.RDS")
fragments<-CreateFragmentObject(path="fragments/fragments.tsv.gz",cells=colnames(intAST))
Fragments(intAST)<-NULL
Fragments(intAST)<-fragments

intMICRO<-readRDS(file="objects/ATAC_microglia.RDS")
fragments<-CreateFragmentObject(path="fragments/fragments.tsv.gz",cells=colnames(intMICRO))
Fragments(intMICRO)<-NULL
Fragments(intMICRO)<-fragments

#read and add fragment path
intEXC<-readRDS(file="objects/ATAC_excitatory.RDS")
fragments<-CreateFragmentObject(path="fragments/fragments.tsv.gz",cells=colnames(intEXC))
Fragments(intEXC)<-NULL
Fragments(intEXC)<-fragments

#read and add fragment path
intL2<-readRDS(file="objects/ATAC_L2-3.RDS")
fragments<-CreateFragmentObject(path="fragments/fragments.tsv.gz",cells=colnames(intL2))
Fragments(intL2)<-NULL
Fragments(intL2)<-fragments

#read and add fragment path
intEXC<-readRDS(file="objects/ATAC_excitatory.RDS")
fragments<-CreateFragmentObject(path="fragments/fragments.tsv.gz",cells=colnames(intEXC))
Fragments(intEXC)<-NULL
Fragments(intEXC)<-fragments

#read and add fragment path
intINH<-readRDS(file="objects/ATAC_inhibitory.RDS")
fragments<-CreateFragmentObject(path="fragments/fragments.tsv.gz",cells=colnames(intINH))
Fragments(intINH)<-NULL
Fragments(intINH)<-fragments


intATAC <- Footprint(
  object = intATAC,
  motif.name = c("CEBPD", "CEBPA", "MAFK"),
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# plot the footprint data for each group of cells
p2 <- PlotFootprint(intATAC, features = c("CEBPD", "CEBPA", "MAFK"))
p2

intOLIGO <- Footprint(
  object = intOLIGO,
  motif.name = c("CEBPD", "CEBPA", "MAFK"),
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# plot the footprint data for each group of cells
p2 <- PlotFootprint(intOLIGO, features = c("CEBPD", "CEBPA", "MAFK"))
p2

intOPC <- Footprint(
  object = intOPC,
  motif.name = c("CEBPD", "CEBPA", "MAFK"),
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# plot the footprint data for each group of cells
p2 <- PlotFootprint(intOPC, features = c("CEBPD", "CEBPA", "MAFK"))
p2

intAST <- Footprint(
  object = intAST,
  motif.name = c("CEBPD", "CEBPA", "MAFK"),
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# plot the footprint data for each group of cells
p2 <- PlotFootprint(intAST, features = c("CEBPD", "CEBPA", "MAFK"))
p2

intEXC <- Footprint(
  object = intEXC,
  motif.name = c("CEBPD", "CEBPA", "MAFK"),
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# plot the footprint data for each group of cells
p2 <- PlotFootprint(intL2, features = c("CEBPD", "CEBPA", "MAFK"))
p2

intL2 <- Footprint(
  object = intL2,
  motif.name = c("CEBPD", "CEBPA", "MAFK"),
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# plot the footprint data for each group of cells
p2 <- PlotFootprint(intL2, features = c("CEBPD", "CEBPA", "MAFK"))
p2

intINH <- Footprint(
  object = intINH,
  motif.name = c("CEBPD", "CEBPA", "MAFK"),
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# plot the footprint data for each group of cells
p2 <- PlotFootprint(intINH, features = c("CEBPD", "CEBPA", "MAFK"))
p2

###################################

#OLIGO
DefaultAssay(intOLIGO)<-"chromvar"
intOLIGO$diagnosis<-intOLIGO$group
Idents(intOLIGO)<-"diagnosis"

oligo.motif.ctrlvC9ALSFTLD <- FindMarkers(
  object = intOLIGO,
  ident.1 = 'C9ALSFTLD',
  ident.2 = 'control',
  only.pos = FALSE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)
save(oligo.motif.ctrlvC9ALSFTLD,file="tables/oligo.motif.ctrlvC9ALSFTLD.rda")
write.csv(oligo.motif.ctrlvC9ALSFTLD,file="tables/oligo.motif.ctrlvC9ALSFTLD.csv")

oligo_motif<-MotifPlot(
  object = intOLIGO,
  motifs = head(rownames(oligo.motif.ctrlvC9ALSFTLD),n=10),
  assay = 'peaks'
)
ggsave('figures/oligo_motif_C9vsCTRL.png',plot=oligo_motif,width=8, height=6, units='in')

oligo_chromvar_C9vsCTRL<-dittoHeatmap(avg_chromvar_oligo,genes=rownames(oligo.motif.ctrlvC9ALSFTLD)[1:100],main = "Oligo DE motifs in C9ALSFTLD",
             annot.by="diagnosis",annot.colors = group_colors, show_colnames = F, show_rownames = F)
ggsave('figures/oligo_chromvar_C9vsCTRL.png', plot=oligo_chromvar_C9vsCTRL, width = 5, height = 6, units='in')


#OPC
DefaultAssay(intOPC)<-"chromvar"
intOPC$diagnosis<-intOPC$group
Idents(intOPC)<-"diagnosis"

OPC.motif.ctrlvC9ALSFTLD <- FindMarkers(
  object = intOPC,
  ident.1 = 'C9ALSFTLD',
  ident.2 = 'control',
  only.pos = F,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

save(OPC.motif.ctrlvC9ALSFTLD,file="tables/OPC.motif.ctrlvC9ALSFTLD.rda")
write.csv(OPC.motif.ctrlvC9ALSFTLD,file="tables/OPC.motif.ctrlvC9ALSFTLD.csv")

OPC_motif<-MotifPlot(
  object = intOPC,
  motifs = head(rownames(OPC.motif.ctrlvC9ALSFTLD),n=10),
  assay = 'peaks'
)
ggsave('figures/OPC_motif_C9vsCTRL.png',plot=OPC_motif,width=8, height=6, units='in')

opc_chromvar_C9vsCTRL<-dittoHeatmap(avg_chromvar_OPC,genes=rownames(diff.act.OPC)[1:100],main = "OPC DE motifs in C9ALSFTLD",
             annot.by="diagnosis",annot.colors = group_colors, show_colnames = F, show_rownames = F)
ggsave('figures/opc_chromvar_C9vsCTRL.png', plot=opc_chromvar_C9vsCTRL, width = 5, height = 6, units='in')

#astrocytes
DefaultAssay(intAST)<-"chromvar"
intAST$diagnosis<-intAST$group
Idents(intAST)<-"diagnosis"

astro.motif.ctrlvC9ALSFTLD <- FindMarkers(
  object = intAST,
  ident.1 = 'C9ALSFTLD',
  ident.2 = 'control',
  only.pos = F,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

save(astro.motif.ctrlvC9ALSFTLD,file="tables/astro.motif.ctrlvC9ALSFTLD.rda")
write.csv(astro.motif.ctrlvC9ALSFTLD,file="tables/astro.motif.ctrlvC9ALSFTLD.csv")

astro_motif<-MotifPlot(
  object = intAST,
  motifs = head(rownames(astro.motif.ctrlvC9ALSFTLD),n=10),
  assay = 'peaks'
)
ggsave('figures/astro_motif_C9vsCTRL.png',plot=astro_motif,width=8, height=6, units='in')

astro_chromvar_C9vsCTRL<-dittoHeatmap(avg_chromvar_astro,genes=rownames(astro.motif.ctrlvC9ALSFTLD)[1:100],main = "Astrocyte DE motifs in C9ALSFTLD",
             annot.by="diagnosis",annot.colors = group_colors, show_colnames = F, show_rownames = F,complex=F)
ggsave('figures/astro_chromvar_C9vsCTRL.png', plot=astro_chromvar_C9vsCTRL, width = 5, height = 6, units='in')


#microglia
DefaultAssay(intMICRO)<-"chromvar"
intMICRO$diagnosis<-intMICRO$group
Idents(intMICRO)<-"diagnosis"

micro.motif.ctrlvC9ALSFTLD <- FindMarkers(
  object = intMICRO,
  ident.1 = 'C9ALSFTLD',
  ident.2 = 'control',
  only.pos = F,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

save(micro.motif.ctrlvC9ALSFTLD,file="tables/micro.motif.ctrlvC9ALSFTLD.rda")
write.csv(micro.motif.ctrlvC9ALSFTLD,file="tables/micro.motif.ctrlvC9ALSFTLD.csv")

micro_motif<-MotifPlot(
  object = intMICRO,
  motifs = head(rownames(micro.motif.ctrlvC9ALSFTLD),n=10),
  assay = 'peaks'
)
ggsave('figures/micro_motif_C9vsCTRL.png',plot=micro_motif,width=8, height=6, units='in')

micro_chromvar_C9vsCTRL<-dittoHeatmap(avg_chromvar_micro,genes=rownames(micro.motif.ctrlvC9ALSFTLD)[1:100], main = "Microglia DE motifs in C9ALSFTLD",
             annot.by="diagnosis",annot.colors = group_colors, show_colnames = F, show_rownames = F)
ggsave('figures/micro_chromvar_C9vsCTRL.png', plot=micro_chromvar_C9vsCTRL, width = 5, height = 6, units='in')

group_colors<-c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2")
CTCF_motif_micro<-VlnPlot_scCustom(intMICRO,features="MA0139.1",group.by = "group",pt.size=0)+
geom_boxplot(width=0.2,fill="white")+
scale_fill_manual(values=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2"))+
ylab("chromvar deviations")+xlab("")+ggtitle("CTCF TF Motif score")+NoLegend()+
theme(plot.title=element_text(face="bold"),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank(),
axis.text.x = element_text(colour="black", size=9), axis.text.y = element_text(colour="black", size=9),
axis.title.x = element_text(colour="black", size=10, face="bold"), axis.title.y = element_text(colour="black", size=10, face="bold"))
ggsave('figures/CTCF_motif_microglia.png', plot=CTCF_motif_micro, width = 5, height = 5, dpi=300, units='in')

Idents(intMICRO)<-"diagnosis"
intMICRO <- Footprint(
  object = intMICRO,
  motif.name = c("CTCF"),
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# plot the footprint data for each group of cells
p2 <- PlotFootprint(intMICRO, features = c("CTCF"),idents=c("control","C9ALSFTLD"))
p2

# excitatory
DefaultAssay(intEXC)<-'chromvar'
Idents(intEXC)<-"group"

exc.motif.ctrlvC9ALSFTLD <- FindMarkers(
  object = intEXC,
  ident.1 = 'C9ALSFTLD',
  ident.2 = 'control',
  only.pos = F,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

save(exc.motif.ctrlvC9ALSFTLD,file="tables/exc.motif.ctrlvC9ALSFTLD.rda")
write.csv(exc.motif.ctrlvC9ALSFTLD,file="tables/exc.motif.ctrlvC9ALSFTLD.csv")

exc_motif<-MotifPlot(
  object = intEXC,
  motifs = head(rownames(exc.motif.ctrlvC9ALSFTLD),n=10),
  assay = 'peaks'
)
ggsave('figures/exc_motif_C9vsCTRL.png',plot=exc_motif,width=8, height=6, units='in')

exc_chromvar_C9vsCTRL<-dittoHeatmap(avg_chromvar_excitatory,genes=rownames(exc.motif.ctrlvC9ALSFTLD)[1:100],main = "Excitatory DE motifs in C9ALSFTLD",
             annot.by="diagnosis",annot.colors = group_colors, show_colnames = F, show_rownames = F)
ggsave('figures/exc_chromvar_C9vsCTRL.png', plot=exc_chromvar_C9vsCTRL, width = 5, height = 6, units='in')


# inhibitory
Idents(intINH)<-"group"
DefaultAssay(intINH)<-'chromvar'

inh.motif.ctrlvC9ALSFTLD <- FindMarkers(
  object = intINH,
  ident.1 = 'C9ALSFTLD',
  ident.2 = 'control',
  only.pos = F,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)
save(inh.motif.ctrlvC9ALSFTLD,file="tables/inh.motif.ctrlvC9ALSFTLD.rda")
write.csv(inh.motif.ctrlvC9ALSFTLD,file="tables/inh.motif.ctrlvC9ALSFTLD.csv")

inh_motif<-MotifPlot(
  object = intINH,
  motifs = head(rownames(inh.motif.ctrlvC9ALSFTLD),n=10),
  assay = 'peaks'
)
ggsave('figures/inh_motif_C9vsCTRL.png',plot=inh_motif,width=8, height=6, units='in')

inh_chromvar_C9vsCTRL<-dittoHeatmap(avg_chromvar_inhibitory,genes=rownames(inh.motif.ctrlvC9ALSFTLD)[1:100],main = "Inhibitory DE motifs in C9ALSFTLD",
             annot.by="diagnosis",annot.colors = group_colors, show_colnames = F, show_rownames = F)
ggsave('figures/inh_chromvar_C9vsCTRL.png', plot=inh_chromvar_C9vsCTRL, width = 5, height = 6, units='in')


inh_chromvar_C9vsCTRL<-dittoHeatmap(intINH,genes=rownames(inh.motif.ctrlvC9ALSFTLD)[1:100],main = "Inhibitory DE motifs in C9ALSFTLD",
             annot.by="diagnosis",annot.colors = group_colors, show_colnames = F, show_rownames = F)
ggsave('figures/inh_chromvar_C9vsCTRL.png', plot=inh_chromvar_C9vsCTRL, width = 5, height = 6, units='in')
