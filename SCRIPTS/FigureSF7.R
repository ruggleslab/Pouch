library(Seurat)
library(ggplot2)
library(matrixStats)
library(gridExtra)
library(RColorBrewer)
library(ggsci)
library(gplots)
library(ComplexHeatmap)
library(circlize)
library(matrixStats)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

colorer <- c("Normal_pouch"="dodgerblue2","Pouchitis"="red3", "UC_colon"="forestgreen", "UC_inflamed"="darkorange1")
colors_clusters = c(pal_d3("category10")(10), pal_d3("category20b")(20), pal_igv("default")(51))


s_obj <- readRDS("OBJECTS/seurat_obj.rds")


gwas_genes <- as.character(read.table("additional_data/IBD_GWAS_genes.txt", F, '\t')[,1])
good_good <- unique(gwas_genes)

good_good <- intersect(good_good, rownames(s_obj))
print(good_good)

counts <- data.frame(s_obj@assays$integrated@scale.data)[good_good,]
tcounts = t(counts)

#
#
#
## 1) What genes are specific for ct
counts <- data.frame(s_obj@assays$integrated@scale.data)[good_good,]
tcounts = t(counts)
agg_counts <- aggregate(tcounts, by=list(s_obj@meta.data$MinorPopulations), "mean")

splitter=as.character(agg_counts$Group.1)

rownames(agg_counts) <- paste0(agg_counts[,1])
agg_counts=agg_counts[,2:ncol(agg_counts)]

de_res <- data.frame(gene=NA, comp = NA, FC = NA)
for(i in 1:length(splitter)){
  x = melt(agg_counts[grep(splitter[i], rownames(agg_counts)),])
  y = melt(agg_counts[grep(splitter[i], rownames(agg_counts), invert=T),])
  y = aggregate(y, by=list(y$variable), "mean")
  
  x$value <- x$value+100
  y$value <- y$value+100
  
  lfc <- (x$value-y$value)/y$value
  
  comp = rep(paste0(splitter[i], " vs all"), nrow(x))
  adder <- data.frame(gene=x$variable, comp = comp, FC = lfc)
  de_res <- rbind(de_res, adder)
}

de_res <- de_res[-1,]
de_res=de_res[order(abs(de_res$FC), decreasing=T),]

namer <- paste0("OBJECTS/GWAS_genes_DE_mean2.txt")
write.table(de_res, namer, row.names=F, sep='\t', quote=F)
#
#
#
#plot it

namer <- paste0("FIGURES/SF7/FigureSF7a_GWAS_genes_rev.pdf")
pdf(namer, height = 8, width = 15)
draw(
  Heatmap(data.matrix(agg_counts), show_row_names = T, show_column_names = T,
          #top_annotation = ha,
          heatmap_legend_param = list(title = "Scaled Value"),
          cluster_rows = T, cluster_columns = T, #row_names_side = 'left',
          column_names_gp = gpar(fontsize=12),
          row_names_gp = gpar(fontsize=12),
          column_names_rot = 90,
          row_title_gp = gpar(fontsize = 10),
          row_names_max_width = unit(10,'cm'),
          use_raster = T,
          cluster_row_slices=F,
          #column_split = splitter,
          col = colorRamp2(c(0, 3), c("white", "red3"))),
  padding = unit(c(2, 2, 3, 2), "cm")
)
dev.off()
