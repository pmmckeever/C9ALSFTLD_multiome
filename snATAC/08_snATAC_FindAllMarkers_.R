library(Seurat)
library(Signac)
library(dplyr)
library(tibble)
library(xlsx)
set.seed(1234)

intATAC<-readRDS(file="objects/intATAC_final.RDS")
fragments<-CreateFragmentObject(path="fragments/fragments.tsv.gz",cells=colnames(intATAC))
Fragments(intATAC)<-NULL
Fragments(intATAC)<-fragments
Idents(intATAC)<-'cellsubtype'


## identify differentially accessible chromatin regions between cellsubtypes - performed similarly to Long et al. Cell Discovery, 2022
DefaultAssay(intATAC) <- "peaks"
Idents(intATAC) <- "cellsubtype"
idents <- as.character(levels(intATAC))
cellsubtype.DARs <- FindAllMarkers(intATAC, 
                                test.use = 'LR',
                                logfc.threshold=0, 
                                min.pct = 0.05,
                                latent.vars = "peak_region_fragments")
cf <- ClosestFeature(intATAC, regions = rownames(cellsubtype.DARs)) # Find the closest feature to a given set of genomic regions
cellsubtype.DARs <- cbind(cellsubtype.DARs, gene=cf$gene_name, gene_biotype = cf$gene_biotype, type = cf$type, distance=cf$distance)
colnames(cellsubtype.DARs)[6:7] <- c("cellsubtype", "genomicRegion")
pi <- Combined.P.FC(cellsubtype.DARs[,c("avg_log2FC", "p_val_adj")], log10P = F)
cellsubtype.DARs$pi <- pi$pi
saveFormat <- lapply(idents, function(x){
  index <- which(cellsubtype.DARs$cellsubtype == x)
  DARs <- cellsubtype.DARs[index,]
  DARs.up <- DARs %>% filter(avg_log2FC>0) %>% arrange(desc(pi))
  DARs.down <- DARs %>% filter(avg_log2FC<0) %>% arrange(pi)
  DARs <- rbind(DARs.up, DARs.down)
  return(DARs)
})
write.xlsx(saveFormat, file = "cellsubtype.all.DARs.xlsx", sheetName = idents, rowNames = F)
saveRDS(cellsubtype.DARs, file = "cellsubtype.all.DARs.rds")

#set a logfc.threshold = 0.25 & p_val_adj < 0.05
cellsubtype.sig.pos.DARs <- cellsubtype.DARs %>% filter(avg_log2FC >=0.25 & p_val_adj < 0.05) %>% arrange(desc(pi)) # 31925 peaks
saveFormat <- lapply(idents, function(x){
  index <- which(cellsubtype.sig.pos.DARs$cellsubtype == x)
  DARs <- cellsubtype.sig.pos.DARs[index,]
  DARs <- DARs %>% arrange(desc(pi))
  return(DARs)
})
names(saveFormat) <- idents
write.xlsx(saveFormat, file = "cellsubtype.sig.DARs.xlsx", sheetName = idents, rowNames = F)
saveRDS(cellsubtype.sig.DARs, file = "cellsubtype.sig.DARs.rds")
