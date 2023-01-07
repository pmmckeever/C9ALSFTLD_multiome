##snRNA/snATAC correlation figures - based on Morabito et al. Nature Genetics 2021 (https://github.com/swaruplab/Single-nuclei-epigenomic-and-transcriptomic-landscape-in-Alzheimer-disease/)

library(Seurat)
library(Signac)
library(reshape2)
library(patchwork)
library(ggpubr)
library(RColorBrewer)
options(stringsAsFactors=FALSE)

intRNA<-readRDS(file="objects/intRNA_final.RDS")
intATAC<-readRDS(file="objects/intATAC_final.RDS")

DefaultAssay(intRNA)<-"RNA"
DefaultAssay(intATAC)<-"RNA"

Idents(intRNA)<-"celltype"
Idents(intATAC)<-"celltype"
celltype<-c('Oligodendrocytes','OPC','Astrocytes','Microglia','Excitatory','Inhibitory')

exp_mat <- GetAssayData(intRNA, slot='data', assay='RNA')
acc_mat <- GetAssayData(intATAC, slot='data', assay='RNA')

genes.use <- intersect(rownames(acc_mat), rownames(exp_mat))

exp_mat <- exp_mat[genes.use,]
acc_mat <- acc_mat[genes.use,]
all.equal(rownames(exp_mat), rownames(acc_mat))
df <- data.frame()

for(cur_celltype in celltype){
  cur_exp_mat <- exp_mat[,intRNA$celltype == cur_celltype]
  cur_acc_mat <- acc_mat[,intATAC$celltype == cur_celltype]
  plot_df <- data.frame(
    atac = rowSums(cur_acc_mat) / ncol(cur_acc_mat),
    rna = rowSums(cur_exp_mat) / ncol(cur_exp_mat),
    gene_name = rownames(exp_mat),
    diagnosis = cur_celltype
  )
  df <- rbind(df, plot_df)
}
save(df, file='data/correlation_RNA_ATAC_celltype.rda')

library(tidyverse)
library(ggrastr)
library(ggpubr)
load('correlation_RNA_ATAC_celltype.rda')
# correlate
p <- ggscatter(
  df, x='rna', y='atac', color='diagnosis', facet.by='diagnosis',
    add = "reg.line",  # Add regression line
 add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
 conf.int = TRUE, # Add confidence interval
 cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
 cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)#+ scale_color_manual(values=unlist(color_scheme_snATAC_celltype)[1:7])
p <- ggplot(df, aes(x=rna, y=atac, color=diagnosis)) +
  rasterise(geom_point(), dpi = 400) +
  geom_smooth(method = 'lm', color='black') +
  stat_cor(p.accuracy = 0.000001, r.accuracy = 0.01, color='black', label.y = 0.9, label.x=0) +
  scale_color_manual(values=unlist(correlation_celltype_colors)[1:7]) +
  ylab('Average Gene Activity (snATAC-seq)') + xlab('Average Gene Expression (snRNA-seq)') + labs_pubr()
pdf('correlate_exp_acc.pdf', width=10, height=7)
p + facet_wrap(~diagnosis, ncol=3)  + theme_pubr()
dev.off()

Idents(intRNA)<-"cellsubtype"
Idents(intATAC)<-"cellsubtype"
cellsubtype<-c('Oligodendrocytes','OPC','AST-FB','AST-PP','Microglia',
  'L2/3','L4','L5/6','L5/6-CC','IN-PV','IN-SST','IN-SV2C','IN-VIP')

exp_mat <- GetAssayData(intRNA, slot='data', assay='RNA')
acc_mat <- GetAssayData(intATAC, slot='data', assay='RNA')

genes.use <- intersect(rownames(acc_mat), rownames(exp_mat))

exp_mat <- exp_mat[genes.use,]
acc_mat <- acc_mat[genes.use,]
all.equal(rownames(exp_mat), rownames(acc_mat))
df <- data.frame()

for(cur_celltype in cellsubtype){
  cur_exp_mat <- exp_mat[,intRNA$cellsubtype == cur_celltype]
  cur_acc_mat <- acc_mat[,intATAC$cellsubtype == cur_celltype]
  # get DEGs for annotation purposes
  #cur_degs <- celltype.markers %>% subset(cluster == cur_celltype & gene %in% genes.use) %>% top_n(5, wt=avg_logFC) %>% .$gene
  plot_df <- data.frame(
    atac = rowSums(cur_acc_mat) / ncol(cur_acc_mat),
    rna = rowSums(cur_exp_mat) / ncol(cur_exp_mat),
    gene_name = rownames(exp_mat),
    diagnosis = cur_celltype
  )
  #plot_df$anno <- ifelse(plot_df$gene %in% cur_degs, plot_df$gene, NA)
  df <- rbind(df, plot_df)
}
save(df, file='data/correlation_RNA_ATAC_cellsubtype.rda')
```

```{r eval=FALSE}
# conda activate r-env
library(tidyverse)
library(ggrastr)
library(ggpubr)
load('correlation_RNA_ATAC_cellsubtype.rda')
#load('color_scheme.rda')
# correlate
p <- ggscatter(
  df, x='rna', y='atac', color='diagnosis', facet.by='diagnosis',
  #label='anno',label.rectangle=TRUE, repel=TRUE,
  add = "reg.line",  # Add regression line
 add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
 conf.int = TRUE, # Add confidence interval
 cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
 cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n")
)#+ scale_color_manual(values=unlist(color_scheme_snATAC_celltype)[1:7])
p <- ggplot(df, aes(x=rna, y=atac, color=diagnosis)) +
  rasterise(geom_point(), dpi = 400) +
  geom_smooth(method = 'lm', color='black') +
  stat_cor(p.accuracy = 0.000001, r.accuracy = 0.01, color='black', label.y = 0.9, label.x=0) +
  scale_color_manual(values=unlist(ATAC_colors_use)[1:13]) +
  ylab('Average Gene Activity (snATAC-seq)') + xlab('Average Gene Expression (snRNA-seq)') + labs_pubr()
pdf('correlate_exp_acc.pdf', width=12, height=6)
p + facet_wrap(~diagnosis, ncol=4)  + theme_pubr()
dev.off()
```

#celltype

```{r eval=FALSE}
library(GeneOverlap)
plot_list <- list()
for(cur_celltype in celltype){
  cur_df <- df %>% subset(diagnosis == cur_celltype)
  cur_genes_atac <- cur_df %>% subset(atac >= quantile(cur_df$atac,0.80)) %>% .$gene
  cur_genes_rna <- cur_df %>% subset(rna <= quantile(cur_df$rna,0.20)) %>% .$gene
  cur_df$overlap <- ifelse(cur_df$gene_name %in% intersect(cur_genes_rna, cur_genes_atac), "Yes", "No")
  data <- as.data.frame(table(cur_df$overlap))
  colnames(data) <- c('category', 'count')
  data$fraction <- data$count / sum(data$count)
  data$ymax <- cumsum(data$fraction)
  data$ymin <- c(0, head(data$ymax, n=-1))
  data$labelPosition <- (data$ymax + data$ymin) / 2
  data$label <- paste0(signif(data$count / sum(data$count) *100, 3), '%')
  # Make the plot
  plot_list[[cur_celltype]] <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect() +
    geom_text( x=2, aes(y=labelPosition, label=label, color=category), size=5) + # x here controls label position (inner / outer)
    scale_fill_manual(values=c(unlist(correlation_celltype_colors[cur_celltype]), 'gray')) +
    scale_color_manual(values=c(unlist(correlation_celltype_colors[cur_celltype]), 'gray')) +
    coord_polar(theta="y") +
    xlim(c(0, 4)) +
    theme_void() + NoLegend()
}
pdf('figures/donut_celltype.pdf', width=10, height=6)
wrap_plots(plot_list, ncol=3)
dev.off()

### cellsubtype
library(GeneOverlap)
plot_list <- list()
for(cur_celltype in cellsubtype){
  cur_df <- df %>% subset(diagnosis == cur_celltype)
  cur_genes_atac <- cur_df %>% subset(atac >= quantile(cur_df$atac,0.80)) %>% .$gene
  cur_genes_rna <- cur_df %>% subset(rna <= quantile(cur_df$rna,0.20)) %>% .$gene
  cur_df$overlap <- ifelse(cur_df$gene_name %in% intersect(cur_genes_rna, cur_genes_atac), "Yes", "No")
  data <- as.data.frame(table(cur_df$overlap))
  colnames(data) <- c('category', 'count')
  data$fraction <- data$count / sum(data$count)
  data$ymax <- cumsum(data$fraction)
  data$ymin <- c(0, head(data$ymax, n=-1))
  data$labelPosition <- (data$ymax + data$ymin) / 2
  data$label <- paste0(signif(data$count / sum(data$count) *100, 3), '%')
  # Make the plot
  plot_list[[cur_celltype]] <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect() +
    geom_text( x=2, aes(y=labelPosition, label=label, color=category), size=5) + # x here controls label position (inner / outer)
    coord_polar(theta="y") +
    xlim(c(0, 4)) +
    theme_void() + NoLegend()
}
pdf('figures/donut_cellsubtype.pdf', width=10, height=6)
wrap_plots(plot_list, ncol=4)
dev.off()
