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
library(reshape)

args = commandArgs(trailingOnly=TRUE)

colorer <- c("Normal_pouch"="dodgerblue2","Pouchitis"="red3", "UC_colon"="forestgreen", "UC_inflamed"="darkorange1")
colors_clusters = c(pal_d3("category10")(10), pal_d3("category20b")(20), pal_igv("default")(51))

big_colorer <- c(cycling_b_cells=colors_clusters[13], 
                 follicular_cells=colors_clusters[14], 
                 gc_cells=colors_clusters[15], 
                 plasma_NFKBIA_hi=colors_clusters[16], 
                 plasma_NFKBIA_lo=colors_clusters[17])

level_fix = c("cycling_b_cells", 
              "follicular_cells", 
              "gc_cells", 
              "plasma_NFKBIA_hi", 
              "plasma_NFKBIA_lo")

s_obj = readRDS("OBJECTS/B_cells/seurat_obj.rds")

s_obj@meta.data$MinorPopulations = factor(s_obj@meta.data$MinorPopulations, levels = level_fix)
Idents(s_obj) = s_obj@meta.data$MinorPopulations

colorer <- colors_clusters
plot_umap =
  DimPlot(s_obj, reduction = "umap", cells = sample(colnames(s_obj)), pt.size=1, cols.use=big_colorer) +
  theme(aspect.ratio = 1) + scale_color_manual(values=big_colorer)

ggsave(paste0("FIGURES/SF3/FigureSF3a.umap.png"), plot = plot_umap, width = 12, height = 10, units = "in")
ggsave(paste0("FIGURES/SF3/FigureSF3a.umap.pdf"), plot = plot_umap, width = 12, height = 10, units = "in")



#
#
#
# and heatmap of top genes
#
#
#
#
marker_list = read.table("OBJECTS/B_cells/clusters-MinorPopulations-clust5/markers-global/markers.clust5.wilcox.all.csv", T, ',')
clusts = levels(s_obj@meta.data$MinorPopulations)

good_genes <- data.frame(gene=NA, ct = NA)
for(i in 1:length(clusts)){
  type = clusts[i]
  gener = subset(marker_list, cluster == clusts[i])
  gener <- gener[order(gener$avg_logFC, gener$p_val_adj, decreasing=T),]
  
  
  genes = as.character(gener$gene[1:10])
  genes2 = data.frame(gene=genes, ct = clusts[i])
  print(length(genes))
  good_genes <- rbind(good_genes, genes2)
}
#good_genes <- c(good_genes, "CD8A", "CD4", "CD3D", "CD79A", "FOS", "FOXP3", "IL17A")
good_genes <- good_genes[-1,]
good_good <- unique(good_genes$gene)

good_good <- intersect(good_good, rownames(s_obj))
print(good_good)

good_genes2 <- subset(good_genes, good_genes$gene %in% good_good)
good_genes2 <- subset(good_genes2, !is.na(gene))
good_genes2$ct <- factor(good_genes2$ct, levels = levels(s_obj@meta.data$MinorPopulations))
good_genes2 <- good_genes2[order(good_genes2$ct),]

#tab2write
figSF5_tab = subset(marker_list, gene %in% good_good)
write.table(figSF5_tab, "FIGURES/SF5/FigureSF5_Top_markers.txt", sep='\t', row.names=F, quote=F)
#

counts <- data.frame(s_obj@assays$integrated@scale.data)[good_genes2$gene,]
tcounts = t(counts)

rownames(s_obj@meta.data) <- make.names(colnames(s_obj))
clusterofchoice = "MinorPopulations"
cluster2 = s_obj@meta.data[colnames(counts),clusterofchoice]
cluster2 = factor(cluster2, levels = levels(s_obj@meta.data$MinorPopulations))

#scounts <- aggregate(tcounts, by=list(cluster2), 'median')
#rownames(scounts) <- scounts[,1]
#scounts <- scounts[,-1]

#scounts[is.na(scounts)] <- 0

ha <- rowAnnotation(df = data.frame(Cluster=cluster2),
                    col=list(Cluster = big_colorer))

namer <- paste0("FIGURES/SF3/FigureS3B_heatmap.png")
png(namer, height = 10, width = 15, units = "in", res=150)
Heatmap(data.matrix(tcounts), show_row_names = F, show_column_names = T,
        #top_annotation = ha,
        heatmap_legend_param = list(title = "Scaled Value"),
        cluster_rows = T, cluster_columns = F, #row_names_side = 'left',
        column_names_gp = gpar(fontsize=12),
        row_title_gp = gpar(fontsize = 10),
        row_names_max_width = unit(10,'cm'),
        use_raster = T,
        cluster_row_slices=F,
        split = cluster2,
        left_annotation = ha,
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
dev.off()


#
#
#
#
##### bxoplots and stats

colors_samples = c(brewer.pal(5, "Set1"), brewer.pal(8, "Dark2"), pal_igv("default")(51))
colors_clusters = c(pal_d3("category10")(10), pal_d3("category20b")(20), pal_igv("default")(51))

colorer <- c("Normal_pouch"="dodgerblue2","Pouchitis"="red3", "UC_colon"="forestgreen", "UC_inflamed"="darkorange1")

cluster = s_obj@meta.data[,"MinorPopulations"]
Idents(s_obj) <- s_obj@meta.data$pouch.status

barplot = data.frame(PouchStatus=Idents(s_obj), Cluster = cluster)
barplot$Sample <- s_obj@meta.data$orig.ident
bb=ggplot(barplot, aes(PouchStatus, fill=cluster)) + geom_bar() + scale_fill_manual(values=colors_clusters) + theme_bw() + ylab("Cell count")
bb2=ggplot(barplot, aes(PouchStatus, fill=cluster)) + geom_bar(position="fill") + scale_fill_manual(values=colors_clusters) + theme_bw() + ylab("Proportion")

df = table(barplot$Cluster, barplot$Sample, barplot$PouchStatus)
dfm = melt(df)
dfm = subset(dfm, value > 0)
colnames(dfm) <- c("Cluster", "Sample", "PouchStatus", "value")
dfm$Cluster <- factor(dfm$Cluster, levels = level_fix)
dfm$PouchStatus <- factor(dfm$PouchStatus, levels = c("Normal_pouch", "Pouchitis", "UC_inflamed"))
melt_plot1 <- ggplot(dfm, aes(PouchStatus, value, color=PouchStatus)) + 
  geom_boxplot(outlier.shape=NA) + geom_jitter(width = 0.2) + xlab("PouchStatus") +
  theme_bw() + scale_color_manual(values = colorer) + ylab("Cell Count") + 
  theme(axis.text.x=element_blank()) +
  facet_wrap(~Cluster, nrow=1, scales="free_y")

dfm$totals <- 0
samps <- as.character(unique(dfm$Sample))
for(i in 1:length(samps)){
  curr = subset(s_obj@meta.data, orig.ident == samps[i])
  totaler = as.numeric(curr$cell_total[1])
  dfm$totals[dfm$Sample == samps[i]] <- totaler
}

dfm$totals <- as.numeric(dfm$totals)
dfm$perc <- (dfm$value/dfm$totals)*100
write_to <- paste0("OBJECTS/B_cells/clusters-resolutions/clust_box_MinorPopulations_values.txt")
write.table(dfm, write_to, sep='\t', row.names=F, quote=F)

dfm$Cluster <- factor(dfm$Cluster, levels = level_fix)

melt_plot2 <- ggplot(dfm, aes(PouchStatus, perc, color=PouchStatus)) + 
  geom_boxplot(outlier.shape=NA) + geom_jitter(width = 0.2) + xlab("PouchStatus") +
  theme_bw() + scale_color_manual(values = colorer) + ylab("Cell Percentage") + 
  theme(axis.text.x=element_blank()) +
  facet_wrap(~Cluster, nrow=1, scales="free_y")

bb3 = ggplot(barplot, aes(cluster, fill=PouchStatus)) + geom_bar() + scale_fill_manual(values=colorer) + theme_bw() + ylab("Cell Count")
bb4 = ggplot(barplot, aes(cluster, fill=PouchStatus)) + geom_bar(position="fill") + scale_fill_manual(values=colorer) + theme_bw() + ylab("Proportion")
bb = arrangeGrob(bb,bb2,bb3,bb4, ncol=2)

#write_to <- list.files(args[1])[grep(gsub("\\.","", args[2]), list.files(args[1]))]
write_to <- paste0("OBJECTS/B_cells/clusters-resolutions/clust_bar_MinorPopulations_values.pdf")

ggsave(write_to, plot = bb, width = 10, height = 10, units = "in")

mm = arrangeGrob(melt_plot1, melt_plot2, ncol=1)
write_to <- paste0("FIGURES/SF3/FigureSF3C_MinorPopulations_boxplot.pdf")
ggsave(write_to, plot = mm, width = 10, height = 6, units = "in")



#
#
#
#
##
###### LogFC heatmap between pouch conditions
library(reshape)

de_table = read.table("OBJECTS/B_cells/clusters-MinorPopulations-clust5/diff-expression-orig.ident/de.clust5.orig.ident.wilcox.all.csv", T, sep=',')

cells = as.character(unique(de_table$cluster))
for(i in 1:length(cells)){
  cell_table=subset(de_table, cluster==cells[i])
  cell_table_de=subset(cell_table, abs(avg_logFC)>0.6 & p_val_adj<0.05)
  good_genes = as.character(unique(cell_table_de$gene))
  cell_table_filt = subset(cell_table, gene %in% good_genes)
  
  if(nrow(cell_table_filt)>0){
    de_mat = cast(cell_table_filt, gene ~ group1+group2, value='avg_logFC')
    rownames(de_mat) = de_mat[,1]
    de_mat=data.matrix(de_mat[,-1])
    
    de_mat[,1]=de_mat[,1]*-1
    de_mat[,3]=de_mat[,3]*-1
    colnames(de_mat) = c("Pouchitis vs. Normal Pouch", "UC Inflamed vs. Normal Pouch",
                         "Pouchitis vs. UC Inflamed")
    de_mat[is.na(de_mat)]<-0
    if(nrow(de_mat)>50){row_size=4} else{row_size=8}
    
    namer <- paste0("OBJECTS/B_cells/clusters-MinorPopulations-clust5/diff-expression-orig.ident/",
                    cells[i], "_FC_heatmap.pdf")
    pdf(namer, height = 10, width = 5)
    print(Heatmap(de_mat, show_row_names = T, show_column_names = T,
                  #top_annotation = ha,
                  heatmap_legend_param = list(title = "Fold Change"),
                  cluster_rows = T, cluster_columns = F, #row_names_side = 'left',
                  column_names_gp = gpar(fontsize=12),
                  row_names_gp = gpar(fontsize=row_size),
                  row_title_gp = gpar(fontsize = 10),
                  #row_names_max_width = unit(10,'cm'),
                  use_raster = T,
                  #cluster_row_slices=F,
                  #split = cluster2,
                  #left_annotation = ha,
                  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
    )
    dev.off()
  }
  
}

