library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(GenomeInfoDb)
library(monocle3)
library(cicero)
library(SeuratWrappers)

intATAC<-readRDS(file="objects/ATAC_AST-PP.RDS")

Idents(intATAC)<-"diagnosis"

DefaultAssay(intATAC)<-"ATAC"

fragments<-CreateFragmentObject(path="fragments/fragments.tsv.gz",cells=colnames(intATAC))
Fragments(intATAC)<-NULL
Fragments(intATAC)<-fragments

load('objects/genome_no_internet/genome_ASTPP.rda')

intATAC.cds<-as.cell_data_set(intATAC)
intATAC.cic<-make_cicero_cds(intATAC.cds,reduced_coordinates=reducedDims(intATAC.cds)$UMAP)
conns<-run_cicero(intATAC.cic,genomic_coords=genome.df,sample_num=100)
head(conns)

ccans<-generate_ccans(conns)
head(ccans)

save(conns, file = "objects/cicero/ASTPP_cicero_conns.rda")
save(ccans, file = "objects/cicero/ASTPP_cicero_ccans.rda")

links<-ConnectionsToLinks(conns=conns,ccans=ccans)

saveRDS(links,file="objects/cicero/intATAC_cicero_links_unadded_ASTPP.RDS")

Links(intATAC)<-links

saveRDS(intATAC,file="objects/cicero/intATAC_cicero_links_ASTPP.RDS")

pdf("figures/CoveragePlot_C9orf72_ASTPP.pdf", height=7, width=11)
CoveragePlot(
  object = intATAC,
  region = "chr9-27546545-27573866",
  extend.upstream = 3000,
  extend.downstream = 3000,
  ncol = 1
)
dev.off()